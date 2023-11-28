

var inBands = ["B2","B3","B4","B5","B6","B8","B11"];
var outBands = inBands.concat("NDVI","NDWI","NDBI","BSI");

var lidar = ee.ImageCollection(image).merge(image2).merge(image3).merge(image4).merge(image5).merge(image6).mosaic();
var lidarclip = lidar.gte(1.5)
  // var training1 = table3.filterBounds(table2)


Map.centerObject(roi1,13);


function maskS2clouds(image) {
  var qa = image.select('QA60');
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask).divide(10000);
}
var period_of_interest = ee.Filter.date('2015-7-01', '2015-10-31');
var dataset = ee.ImageCollection('COPERNICUS/S2')
                  .filter(period_of_interest)
                  .filterBounds(roi1)
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',10))
                  .map(maskS2clouds);
var base_bands = dataset.mean().select(["B2","B3","B4","B5","B6","B8"])
var addNDVI = function(image) {var ndvi = image.normalizedDifference(['B8', 'B4'])
  .rename('NDVI')
  .copyProperties(image,['system:time_start']);
  return image.addBands(ndvi);
};

var addBSI = function(image) {var bsi = image.expression(
  '((RED + SWIR) - (NIR + BLUE)) / ((RED + SWIR) + (NIR + BLUE)) ', 
  {
    'RED': image.select('B4'), 
    'BLUE': image.select('B2'),
    'NIR': image.select('B8'),
    'SWIR': image.select('B11'),
  }
)
  .rename('BSI')
  .copyProperties(image,['system:time_start']);
  return image.addBands(bsi);
};


var addNDBI = function(image) {var ndbi = image.normalizedDifference(['B11', 'B8'])
  .rename('NDBI')
  .copyProperties(image,['system:time_start']);
  return image.addBands(ndbi);
};
var addNDRE1 = function(image){var ndre1 = image.expression(
  '((NIR - RE) / (NIR + RE))  ', 
  {
    'NIR': image.select('B8'), 
    
    'RE': image.select('B5'),
  }
).rename('NDRE1')
return image.addBands(ndre1)
}
var addNDRE2 = function(image){var ndre2 = image.expression(
  '((NIR - RE) / (NIR + RE))  ', 
  {
    'NIR': image.select('B8'), 
    
    'RE': image.select('B6'),
  }
).rename('NDRE2')
return image.addBands(ndre2)
}

// var ndre = rapideye.expression(
//   '((NIR - RE) / (NIR + RE))  ', 
//   {
//     'NIR': rapideye.select('b5'), 
    
//     'RE': rapideye.select('b4'),
//   }
// ).rename('NDRE')
// var rededge = rapideye.select('b4').rename('REDEDGE')

// var gcc = rapideye.expression(
//   '((green) / (blue + green+ red))  ', 
//   {
//     'red': rapideye.select('b1'), 
    
//     'green': rapideye.select('b2'),
    
    
//     'blue': rapideye.select('b3'),
//   }
// ).rename('GCC')

var addevi = function(image){var evi = image.expression(
    '2.5 * ((B08 - B04) / ((B08 + 6.0 * B04 - 7.5 * B02) + 1.0))', {
      'B08': image.select('B8'),
      'B04': image.select('B4'),
      'B02': image.select('B2'),
}).rename('EVI')
return image.addBands(evi)
}
var addgccsen = function(image){var gcc_sen = image.expression(
  '((green) / (blue + green+ red))  ', 
  {
    'red': image.select('B4'), 
    
    'green': image.select('B3'),
    
    
    'blue': image.select('B2'),
  }
).rename('GCC_SEN')
return image.addBands(gcc_sen)
}
//Collection with all images also containing the NDVI and BSI indices
var collection = dataset.map(addBSI)
.map(addNDVI).map(addNDBI)
.map(addNDRE1).map(addNDRE2).map(addgccsen).map(addevi);
// print(collection,'Collection with inBands and Statistic indices');

// //min
// Define and calculate the median bands and the other index statistics
var band_mosaic = collection.select(inBands).mosaic();
var EVI = collection.select('EVI').reduce(ee.Reducer.mean()).rename("EVI")

var NDRE1 = collection.select('NDRE1').reduce(ee.Reducer.mean()).rename("NDRE1")
var GCC_SEN = collection.select('GCC_SEN').reduce(ee.Reducer.mean()).rename("GCC_SEN")

var ndvi = collection.select('NDVI').reduce(ee.Reducer.mean()).rename("NDVI");

// var ndwimax = collection.select('NDWI').reduce(ee.Reducer.max()).rename("NDWI_MAX");
// var ndwimean = collection.select('NDWI').reduce(ee.Reducer.mean()).rename("NDWI");
// var ndwistd = collection.select('NDWI').reduce(ee.Reducer.stdDev()).float().rename("NDWI_STD");
// var ndbimax = collection.select('NDBI').reduce(ee.Reducer.max()).rename("NDBI_MAX");
var ndbi = collection.select('NDBI').reduce(ee.Reducer.mean()).rename("NDBI");
var ndbistd = collection.select('NDBI').reduce(ee.Reducer.stdDev()).float().rename("NDVI_STD")
var compclip = band_mosaic.addBands(EVI).addBands(GCC_SEN).addBands(lidar).addBands(ndvi).addBands(ndbi);
// Map.addLayer(compclip.select(['NDBI','NDVI','EVI','GCC','B2','B3','B5','B4']))
var ndbiclip = compclip.select('NDBI').gte(-0.01);
var gccclip = compclip.select('GCC_SEN').lte(0.35);
var ndviclip = compclip.select('NDVI').lte(0.3);
var ndviclip1 = compclip.select('NDVI').lte(0.5);
var har = lidarclip.subtract(ndviclip1)
var ndviclip1 = lidarclip.subtract(ndviclip).subtract(ndbiclip).subtract(gccclip)
var har1 = har.gte(0.5)
var trial = lidar.gt(1.5).and(compclip.select('NDVI_RE').gte(0.4))
.and(compclip.select('NDBI').lte(0.16)).and(compclip.select('GCC').gte(0.39))
var final_tree = ndviclip1.gte(0.5)
// Map.addLayer(final_tree.clip(roi),{},'Tree_crowns')

// Map.addLayer(trial.clip(roi),{},'Tree_crownssssss')

// print("Composition", compclip);
var final_tree_bands = final_tree.addBands(base_bands).addBands(ndvi).addBands(lidarclip)
.addBands(NDRE1).addBands(GCC_SEN).addBands(EVI)
// .addBands(ndvire).addBands(ndre).addBands(rededge).addBands(gcc);
Map.addLayer(final_tree_bands.clip(roi1),{},'Indices_Image',false)
// var ndvimin = final_tree_bands.select('NDVI').reduce(ee.Reducer.min()).rename("NDVI_MIN");
// var composite = final_tree_bands.addBands('NDVI_MIN')
var img1 = final_tree_bands.select('B8','B5','b1','EVI','NDVI','GCC_SEN','NDRE1')
// ,'NDVI_RE','NDRE','REDEDGE','GCC','EVI')
// Map.addLayer(final_bands.clip(table),{})
// var final_bands= compclip.select("NDVI_MIN",'b1');


// var final_tree_bands1 = final_tree.addBands(ndvimin).addBands(ndvimax).addBands(bsimean).addBands(ndvi).addBands(bsimax).addBands(lidarclip)
// .addBands(NDRE1).addBands(NDRE2).addBands(GCC_SEN)
// .addBands(ndvire).addBands(ndre).addBands(rededge).addBands(gcc).addBands(EVI);
// print("Composition", compclip1);

// var img1= final_tree_bands1.select('BSI','BSI_MAX','b1','NDVI_MAX',"NDVI_MIN",'NDVI'
// ,'NDRE1','NDRE2','GCC_SEN').clip(table).divide(255)
// ,'NDVI_RE','NDRE','REDEDGE','GCC','EVI');
// final_bands1.bandNames()
// print(final_bands1)
// var final_bands2= compclip.select("NDVI_MAX");
// Map.addLayer(final_bands2, {}, "Finalmax");



var seeds = ee.Algorithms.Image.Segmentation.seedGrid(36,'hex');

// Run SNIC on the regular square grid.
var snic = ee.Algorithms.Image.Segmentation.SNIC({
  image: img1, 
  size: 11,
  compactness: 11,
  connectivity: 4,
  neighborhoodSize:0,
  seeds: seeds
}).select(['B8_mean','B5_mean','b1_mean','EVI_mean','NDVI_mean','GCC_SEN_mean','NDRE1_mean', 'clusters'],
['B8','B5','b1','EVI','NDVI','GCC_SEN','NDRE1', 'clusters'])
var bandNames = snic.bandNames();
print('Band names:', bandNames);
var clusters = snic.select('clusters')
Map.addLayer(clusters.randomVisualizer(), {}, 'clusters',false)
// Map.addLayer(snic, {bands: ['R', 'G', 'B'], min:0, max:1, gamma: 0.8}, 'means', false)

// // Compute per-cluster stdDev.
var stdDev = img1.addBands(clusters).reduceConnectedComponents(ee.Reducer.stdDev(), 'clusters', 256)
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
var training_bands = ['B8','B5','b1','EVI','NDVI','GCC_SEN','NDRE1', 'clusters']
var objectPropertiesImage = ee.Image.cat([
  snic.select(training_bands),
  stdDev,
  // area,
  perimeter,
  width,
  height
]).float();
 
var training_data = objectPropertiesImage.sampleRegions({
  collection: training2,
  properties: ['Class'],
  
  scale: 1
})
var trainingData = training_data.randomColumn();
var trainset = trainingData.filter(ee.Filter.lessThan('random', 0.8));
var testset = trainingData.filter(ee.Filter.greaterThanOrEquals('random', 0.8));
var landcoverpalette = [
  'FFF300',//coniferous
  'FF4900'];//deciduous
  
  var classifi =  ee.Classifier.libsvm(
  // {
  // decisionProcedure: 'Voting',
  // svmType: 'C_SVC',
  // kernelType: 'POLY',
  // shrinking: false,
  // // degree: 3, //Valid only for 'POLY' kernelType
  // gamma: null, //Valid for 'POLY', 'RBF' and 'SIGMOID' kernelType - null:defoult value= reciprocal of the number of features
  // coef0: 1, //Valid for 'POLY' and 'SIGMOID' kernelType - is the value of coeff0 in kernel function. Default=0
  // cost: 0.5}
  ).train({
  features: trainset,
  classProperty: 'Class',
  inputProperties: objectPropertiesImage.bandNames()
});
  var classifie =  ee.Classifier.smileCart().train({
  features: trainset,
  classProperty: 'Class',
  inputProperties: objectPropertiesImage.bandNames()
});
// ee.Classifier.libsvm({
//     decisionProcedure: 'voting',
//   kernelType: 'RBF',
//   gamma: 1,
//   cost: 10
// }).train(trainset, 'Class',training_bands);
// var classifier_svm = classifi.train(training, 'class', training_bands);
var classifier = ee.Classifier.smileRandomForest(200,2).train(trainset, 'Class', objectPropertiesImage.bandNames())
Map.addLayer(objectPropertiesImage.classify(classifier).clip(table2),{palette: landcoverpalette, min:0, max: 1}, 'Classified objects RF')
var rf_image = objectPropertiesImage.classify(classifier).clip(table2).visualize({
  min: 0,
  max: 1,
  palette: ['FFF300','FF4900']
});
Export.image.toAsset({
  image: rf_image.clip(table2),
  description: 'object_se_rf_2021',
  scale:1,
  region: roi1, 
  maxPixels:10e12
});
Map.addLayer(objectPropertiesImage.classify(classifi).clip(table2),{palette: landcoverpalette, min:0, max: 1}, 'Classified objects SVM')
var svm_image = objectPropertiesImage.classify(classifier).clip(table2).visualize({
  min: 0,
  max: 1,
  palette: ['FFF300','FF4900']
});
Export.image.toAsset({
  image: svm_image.clip(table2),
  description: 'object_se_svm_2021',
  scale:1,
  region: roi1, 
  maxPixels:10e12
});
Map.addLayer(objectPropertiesImage.classify(classifie).clip(table2),{palette: landcoverpalette, min:0, max: 1}, 'Classified objects CART')
var cart_image = objectPropertiesImage.classify(classifier).clip(table2).visualize({
  min: 0,
  max: 1,
  palette: ['FFF300','FF4900']
});
Export.image.toAsset({
  image: cart_image.clip(table2),
  description: 'object_se_cart_2021',
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






