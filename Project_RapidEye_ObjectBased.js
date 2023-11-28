


Map.centerObject(roi1,15)
Map.addLayer(rapideye.clip(roi1),{min:67.0084,max:1993.97,bands: ["b5","b4","b2"]},'Rapid eye',false)
// var training1 = ee.FeatureCollection(Coniferous).merge(Deciduous)

var ndvire = rapideye.expression(
  '((NIR - RED) / (NIR + RED))  ', 
  {
    'NIR': rapideye.select('b5'), 
    
    'RED': rapideye.select('b1'),
  }
).rename('NDVI_RE')
var ndre = rapideye.expression(
  '((NIR - RE) / (NIR + RE))  ', 
  {
    'NIR': rapideye.select('b5'), 
    
    'RE': rapideye.select('b4'),
  }
).rename('NDRE')
var rededge = rapideye.select('b4').rename('REDEDGE')

var gcc = rapideye.expression(
  '((green) / (blue + green+ red))  ', 
  {
    'red': rapideye.select('b1'), 
    
    'green': rapideye.select('b2'),
    
    
    'blue': rapideye.select('b3'),
  }
).rename('GCC')
var evi = rapideye.expression(
    '2.5 * ((B08 - B04) / ((B08 + 6.0 * B04 - 7.5 * B02) + 1.0))', {
      'B08': rapideye.select('b5'),
      'B04': rapideye.select('b1'),
      'B02': rapideye.select('b3'),
}).rename('EVI')
var lidar = ee.ImageCollection(image).merge(image2).merge(image3).merge(image4).merge(image5).merge(image6).mosaic().clip(roi1);
var lidar = lidar.select('b1').rename('b6')

var rapideye_indices = rapideye.addBands(ndvire).addBands(ndre).addBands(gcc).addBands(lidar).addBands(evi)
// .addBands(ndvimin).addBands(bsimean).addBands(ndvi).addBands(EVI)
var bands = ['b1', 'b2', 'b3', 'b4','b5','b6']
var img = ee.Image(rapideye_indices).select(['b1','b2','b3','b4','b5','b6','NDVI_RE','NDRE','GCC','EVI']).clip(table).divide(255)
print(img.bandNames())
Map.centerObject(roi1, 13)
Map.addLayer(img, {gamma: 0.8}, 'RGBN', false)

var seeds = ee.Algorithms.Image.Segmentation.seedGrid(36,'hex');

// Run SNIC on the regular square grid.
var snic = ee.Algorithms.Image.Segmentation.SNIC({
  image: img, 
  size: 11,
  compactness: 11,
  connectivity: 4,
  neighborhoodSize:0,
  seeds: seeds
}).select(['b4_mean','b5_mean','b6_mean','NDVI_RE_mean','NDRE_mean','GCC_mean','EVI_mean', 'clusters'],
['RE','N','elevation','NDVI','NDRE','GCC','EVI', 'clusters'])
var bandNames = snic.bandNames();
print('Band names:', bandNames);
var clusters = snic.select('clusters')
Map.addLayer(clusters.randomVisualizer(), {}, 'clusters',false)
// Map.addLayer(snic, {bands: ['R', 'G', 'B'], min:0, max:1, gamma: 0.8}, 'means', false)

// // Compute per-cluster stdDev.
var stdDev = img.addBands(clusters).reduceConnectedComponents(ee.Reducer.stdDev(), 'clusters', 256)
Map.addLayer(stdDev, {min:0, max:0.1}, 'StdDev', false)

// Area, Perimeter, Width and Height
var area = ee.Image.pixelArea().addBands(clusters).reduceConnectedComponents(ee.Reducer.sum(), 'clusters', 256)
Map.addLayer(area, {min:50000, max: 500000}, 'Cluster Area', false)

var minMax = clusters.reduceNeighborhood(ee.Reducer.minMax(), ee.Kernel.square(1));

var perimeterPixels = minMax.select(0).neq(minMax.select(1)).rename('perimeter');

Map.addLayer(perimeterPixels, {min: 0, max: 1}, 'perimeterPixels',false);

var perimeter = perimeterPixels.addBands(clusters)
    .reduceConnectedComponents(ee.Reducer.sum(), 'clusters', 256);
Map.addLayer(perimeter, {min: 100, max: 400}, 'Perimeter size', false);

var sizes = ee.Image.pixelLonLat().addBands(clusters).reduceConnectedComponents(ee.Reducer.minMax(), 'clusters', 256)
var width = sizes.select('longitude_max').subtract(sizes.select('longitude_min')).rename('width')
var height = sizes.select('latitude_max').subtract(sizes.select('latitude_min')).rename('height')
Map.addLayer(width, {min:0, max:0.02}, 'Cluster width', false)
Map.addLayer(height, {min:0, max:0.02}, 'Cluster height', false)
var training_bands = ['N','RE','NDVI','NDRE','GCC']
var objectPropertiesImage = ee.Image.cat([
  snic,
  stdDev,
  // area,
  perimeter,
  width,
  height
]).float();

var trainingbands = ['N','RE','NDVI','NDRE','GCC','width','height']
var training_data = objectPropertiesImage.sampleRegions({
  collection: table4,
  properties: ['Class'],
  
  scale: 1
})
var trainingData = training_data.randomColumn();
var trainset = trainingData.filter(ee.Filter.lessThan('random', 0.8));
var testset = trainingData.filter(ee.Filter.greaterThanOrEquals('random', 0.8));
var landcoverpalette = [
  'FFF300',//coniferous
  'FF4900'];//deciduous
  
  var classifi =  ee.Classifier.libsvm().train({
  features: trainset,
  classProperty: 'Class',
  inputProperties: trainingbands
});
  var classifie =  ee.Classifier.smileCart().train({
  features: trainset,
  classProperty: 'Class',
  inputProperties: objectPropertiesImage.bandNames()
});

var classifier = ee.Classifier.smileRandomForest(450,2).train(trainset, 'Class', objectPropertiesImage.bandNames())
Map.addLayer(objectPropertiesImage.classify(classifier).clip(table2),{palette: landcoverpalette, min:0, max: 1}, 'Classified objects RF')
var rf_image = objectPropertiesImage.classify(classifier).clip(table2).visualize({
  min: 0,
  max: 1,
  palette: ['FFF300','FF4900']
});
Export.image.toAsset({
  image: rf_image,
  description: 'object_re_rf',
  scale:1,
  region: roi1, 
  maxPixels:10e12
});
Map.addLayer(objectPropertiesImage.classify(classifi).clip(table2),{palette: landcoverpalette, min:0, max: 1}, 'Classified objects SVM')
var svm_image = objectPropertiesImage.classify(classifi).clip(table2).visualize({
  min: 0,
  max: 1,
  palette: ['FFF300','FF4900']
});
Export.image.toAsset({
  image: svm_image,
  description: 'object_re_svm',
  scale:1,
  region: roi1, 
  maxPixels:10e12
});
Map.addLayer(objectPropertiesImage.classify(classifie).clip(table2),{palette: landcoverpalette, min:0, max: 1}, 'Classified objects CART')
var cart_image = objectPropertiesImage.classify(classifie).clip(table2).visualize({
  min: 0,
  max: 1,
  palette: ['FFF300','FF4900']
});
Export.image.toAsset({
  image: cart_image,
  description: 'object_re_cart',
  scale:1,
  region: roi1, 
  maxPixels:10e12
});
var confusionmatrix = ee.ConfusionMatrix(testset.classify(classifier).errorMatrix({actual: 'Class', predicted: 'classification'}));

print('confusion Matrix:RF', confusionmatrix);
print('Overall Accuracy:RF', confusionmatrix.accuracy());
print('Overall kappa:RF', confusionmatrix.kappa());
print("Consumer's accuracy:RF", confusionmatrix.consumersAccuracy());
print("Producer's accuracy:RF", confusionmatrix.producersAccuracy());
var confusionmatrix = ee.ConfusionMatrix(testset.classify(classifi).errorMatrix({actual: 'Class', predicted: 'classification'}));

print('confusion Matrix:SVM', confusionmatrix);
print('Overall Accuracy:SVM', confusionmatrix.accuracy());
print('Overall kappa:SVM', confusionmatrix.kappa());
print("Consumer's accuracySVM", confusionmatrix.consumersAccuracy());
print("Producer's accuracySVM", confusionmatrix.producersAccuracy());
var confusionmatrix = ee.ConfusionMatrix(testset.classify(classifie).errorMatrix({actual: 'Class', predicted: 'classification'}));

print('confusion Matrix:CART', confusionmatrix);
print('Overall Accuracy:CART', confusionmatrix.accuracy());
print('Overall kappa:CART', confusionmatrix.kappa());
print("Consumer's accuracy:CART", confusionmatrix.consumersAccuracy());
print("Producer's accuracy:CART", confusionmatrix.producersAccuracy());

var legend = ui.Panel({
  style: {
    position: 'bottom-left',
    padding: '8px 15px'
  }
});
 
// Create legend title
var legendTitle = ui.Label({
  value: 'Class',
  style: {
    fontWeight: 'bold',
    fontSize: '18px',
    margin: '0 0 4px 0',
    padding: '0'
    }
});
 
// Add the title to the panel
legend.add(legendTitle);
 
// Creates and styles 1 row of the legend.
var makeRow = function(color, name) {
 
      // Create the label that is actually the colored box.
      var colorBox = ui.Label({
        style: {
          backgroundColor: '#' + color,
          // Use padding to give the box height and width.
          padding: '8px',
          margin: '0 0 4px 0'
        }
      });
 
      // Create the label filled with the description text.
      var description = ui.Label({
        value: name,
        style: {margin: '0 0 4px 6px'}
      });
 
      // return the panel
      return ui.Panel({
        widgets: [colorBox, description],
        layout: ui.Panel.Layout.Flow('horizontal')
      });
};
 
//  Palette with the colors
var palette =['FFF300',//coniferous
  'FF4900'];

// name of the legend
var names = ['Coniferous','Deciduous'];
 
// Add color and and names
for (var i = 0; i < 2; i++) {
  legend.add(makeRow(palette[i], names[i]));
  }  
 
// add legend to map (alternatively you can also print the legend to the console)
Map.add(legend);



 var dict_RF = classifier.explain();
var variable_importance_RF = ee.Feature(null, ee.Dictionary(dict_RF).get('importance'));
  var chart_RF =
ui.Chart.feature.byProperty(variable_importance_RF)
.setChartType('ColumnChart')
.setOptions({
title: 'Random Forest Variable Importance',
legend: {position: 'none'},
hAxis: {title: 'Bands'},
vAxis: {title: 'Importance'}
});
 
print(chart_RF);


var dict_RF = classifie.explain();

 
var variable_importance_RF = ee.Feature(null, ee.Dictionary(dict_RF).get('importance'));
  var chart_RF =
ui.Chart.feature.byProperty(variable_importance_RF)
.setChartType('ColumnChart')
.setOptions({
title: 'Cart Variable Importance',
legend: {position: 'none'},
hAxis: {title: 'Bands'},
vAxis: {title: 'Importance'}
});
 
print(chart_RF);
