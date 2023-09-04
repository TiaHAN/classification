import geemap
import ee
# VPN port
geemap.set_proxy(port=)
try:
    ee.Initialize()
except:
    ee.Authenticate()
    ee.Initialize()
Map = geemap.Map()
Map

# 1）input parameters
region = 'users/wanghan02191127/shanxilinyi'
samp = 'users/wanghan02191127/Shaanxi_6classPoints'
Shaanxi_OrchardAndOther_ClassifiedImg = 'users/wanghan02191127/Shaanxi_OrchardAndOther_ClassifiedImg'
Orchard_ClassifiedImg = Shaanxi_OrchardAndOther_ClassifiedImg.eq(2)
Map.addLayer(Shaanxi_OrchardAndOther_ClassifiedImg,{'min': 1, 'max': 2, 'palette': ['f7fff4', 'd6230c']},'OrchardAndOther_ClassifiedImg')
Map.addLayer(Shaanxi_OrchardAndOther_ClassifiedImg.updateMask(Orchard_ClassifiedImg),{},'OrchardAndOther_ClassifiedImg222')


# 2)defined function
def sentinel2toa(img):
    selected_bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B11', 'B12']
    toa_img = img.select(selected_bands).divide(10000).toDouble()
    qa60_band = img.select(['QA60'])
    scl_band = img.select(['SCL'])
    solar_azimuth = img.get('MEAN_SOLAR_AZIMUTH_ANGLE')
    solar_zenith = img.get('MEAN_SOLAR_ZENITH_ANGLE')
    time_start = img.get('system:time_start')
    toa_img = toa_img.addBands(qa60_band).addBands(scl_band) \
                     .set('solar_azimuth', solar_azimuth) \
                     .set('solar_zenith', solar_zenith) \
                     .set('system:time_start', time_start)
    return toa_img

def cloudMask(toa):
    # Compute several indicators of cloudiness and take the minimum of them
    def rescale(img, thresholds):
        return img.subtract(thresholds[0]).divide(thresholds[1] - thresholds[0])
    # Initialize the cloud score image with a value of 1
    score = ee.Image(1)
    # Clouds are reasonably bright
    score = score.min(rescale(toa.select(['B2']), [0.1, 0.5]))
    score = score.min(rescale(toa.select(['B1']), [0.1, 0.3]))
    score = score.min(rescale(toa.select(['B1']).add(toa.select(['B9'])), [0.15, 0.2]))
    score = score.min(rescale(toa.select(['B4']).add(toa.select(['B3'])).add(toa.select('B2')), [0.2, 0.8]))
    # Clouds are moist
    ndmi = toa.normalizedDifference(['B8A', 'B11'])
    score = score.min(rescale(ndmi, [-0.3, 0.3]))
    # However, clouds are not snow
    ndsi = toa.normalizedDifference(['B3', 'B11'])
    score = score.min(rescale(ndsi, [0.9, 0.4]))
    # Define a cloud score threshold
    cloudScoreThreshold = 0.12
    cloud = score.gt(cloudScoreThreshold)
    mask = cloud.eq(0)
    return toa.updateMask(mask)

# Define the function to add indices
def addIndices(image):
    ndvi = image.normalizedDifference(['B8', 'B4'])
    ndmi = image.normalizedDifference(['B12', 'B3'])
    ndri = image.normalizedDifference(['B11', 'B3'])
    A = image.expression('(RE3-RE2)/RE1', {
        'RE1': image.select('B5'),
        'RE2': image.select('B6'),
        'RE3': image.select('B7')
    })
    return image.addBands(ndvi.rename("NDVI")) \
                .addBands(ndmi.rename('NDMI')) \
                .addBands(ndri.rename("NDRI")) \
                .addBands(A.rename("A"))

# 3) sentinel_2
# 通过分析时序曲线的差异，选择以下波段月际时序曲线作为分类影像
bands = ee.List(['B11', 'B2', 'B7', 'NDVI', 'NDMI', 'NDRI','A'])
s2_data = ee.ImageCollection("COPERNICUS/S2_SR")\
            .filterDate('2021-1-1','2022-1-1')\
             .filterBounds(region)\
             .map(sentinel2toa)\
             .map(cloudMask )\
             .map(addIndices)\
             .select(bands)
             
             
# 4) cmposite 12 month image            
month = ee.List.sequence(1,12,1)
year = ee.Number(2021)
def create_composite(i):
    i = ee.Number(i)
    SDay = ee.Date.fromYMD(year, i, 1)
    EDay = SDay.advance(1, 'month')
    Img = s2_data.filterDate(SDay, EDay).mean()
    time = SDay.advance(15, 'day')
    return Img.set('time', time.format('yyyy-MM-dd'))
month_com_Col = ee.ImageCollection(month.map(create_composite))
# print("month_com_Col:",month_com_Col)


# harmonic regression
time_reference = ee.Date(2021-1-1)

def process_image(img):
    tstamp = ee.Date(img.get('time'))
    tdelta=tstamp.difference(time_reference,'year')
    
    img_fitting=img.select()\
                   .addBands(ee.Image(tdelta.multiply(2*3.141592653589793).cos()).rename('cos'))\
                   .addBands(ee.Image(tdelta.multiply(2*3.141592653589793).cos()).rename('sin'))\
                   .addBands(img.select(bands))\
                   .toDouble()
    return img_fitting    
month_com_col_fit = month_com_Col.map(process_image)          
dependent = bands
harmonicIndependents = ee.List(['constant', 'cos', 'sin'])
harmonic = month_com_col_fit.select(harmonicIndependents.cat(dependent))\
                 .reduce(ee.Reducer.linearRegression(harmonicIndependents.length(), dependent.length()))
coefficients = harmonic.select('coefficients').matrixTranspose()\
                       .arrayFlatten([dependent,harmonicIndependents]).clip(region)
# print('harmonic coeff:',coefficients)


#display regression
def mergeBands(image, previous):
    return ee.Image(previous).addBands(image)
# Define the mapping function
def mapFunc(image):
    def innerMap(e):
        elem = ee.String(e)
        bandname = elem.cat(".*")       
        return image.select(harmonicIndependents) \
                    .multiply(coefficients.select(bandname)) \
                    .reduce('sum') \
                    .rename(elem)   
    newImgCol = ee.List(bands).map(innerMap)    
    newImg = ee.Image(newImgCol.iterate(mergeBands, ee.Image([])))
    return newImg.set('system:time_start', image.get('time'))

# Map over the collection and apply the mapping function
month_com_col_fitted = month_com_col_fit.map(mapFunc)

# Print the result
print(month_com_col_fitted.getInfo())


# 高程数据
dataset = ee.Image('USGS/SRTMGL1_003')   
terrain_data = ee.Algorithms.Terrain(dataset)


# 分类输入影像数据集
Classimg_data =ee.Image(month_com_col_fitted.seclect(bands).iterate(mergeBands,ee.Image([])))\
              .addBands(terrain_data)\
              .updateMask(Orchard_ClassifiedImg)



# 5)RF分类
# 先分样本点在验证
print('先分样本点，在验证')
SampleFeas_Random = samp.randomColumn('random',42)
trainingPoints = SampleFeas_Random.filter(ee.Filter.lt('random',0.7))
testingPoints = SampleFeas_Random.filter(ee.Filter.gte('random',0.7))

trainingPoints_value= Classimg_data.sampleRegions(
  collection = trainingPoints,
  properties = ['Type'],
  scale = 10,
  tileScale =16,
  geometries = True
)
testingPoints_value = Classimg_data.sampleRegions(
  collection= testingPoints,
  properties = ['Type'],
  scale = 10,
  tileScale = 16,
  geometries = True
)
# 分类
AllBandName = Classimg_data.bandNames()
RF_classifier = ee.Classifier.smileRandomForest(numberOfTrees=100, seed=42)
TrainedModel1 = RF_classifier.train(trainingPoints_value, 'Type', AllBandName)
RFclassifiedResult = Classimg_data.select(AllBandName).classify(TrainedModel1)
print('1.2RF分类影像',RFclassifiedResult)

ValidatedResult1= testingPoints_value.classify(TrainedModel1)
ValidatedAccuracy1= ValidatedResult1.errorMatrix('Type', 'classification')

# print('1.3RF分类精度:',
#       'ErrorMatrix:',ValidatedAccuracy1,'Validate accuracy:',ValidatedAccuracy1.accuracy(),
#       'Validate accuracy:',ValidatedAccuracy1.accuracy(),'User acc:',ValidatedAccuracy1.consumersAccuracy(),
#       'Prod acc:',ValidatedAccuracy1.producersAccuracy(),'Kappa:',ValidatedAccuracy1.kappa())

# 导出影像
outimg = RFclassifiedResult.toByte().remap([0,1,2,3,4,5],[1,2,3,4,5,6])
task = ee.image.toDrive(
  image = outimg.updateMask(Orchard_ClassifiedImg), 
  description = 'Shaanxi_6ClassImage_RF',
  folder = 'AppleMap', 
  region = region,
  scale = 10,
  maxPixels = 1e13
)

# 生成概率密度结果
classProperty= 'Type'
RF_classifier_Multipro_fun = RF_classifier.setOutputMode('MULTIPROBABILITY')
trained_classifier = RF_classifier_Multipro_fun.train(
    features=trainingPoints_value,
    classProperty = classProperty,
    inputProperties=AllBandName)

RF_Multipro = Classimg_data.classify(RF_classifier_Multipro_fun)
RF_Multipro_img = (ee.Image.cat([
    RF_Multipro.arraySlice(0, 0, 1).arrayProject([0]).arrayFlatten([['Apple_pro']]),
    RF_Multipro.arraySlice(0, 1, 2).arrayProject([0]).arrayFlatten([['Pear_pro']]),
    RF_Multipro.arraySlice(0, 2, 3).arrayProject([0]).arrayFlatten([['Peach_pro']]),
    RF_Multipro.arraySlice(0, 3, 4).arrayProject([0]).arrayFlatten([['Kiwi_pro']]),
    RF_Multipro.arraySlice(0, 4, 5).arrayProject([0]).arrayFlatten([['Grape_pro']]),
    RF_Multipro.arraySlice(0, 5, 6).arrayProject([0]).arrayFlatten([['Apricot_pro']])
]))
print("1.4RF概率密度影像:",RF_Multipro_img)


# 6）k_means
trainingFea_Kmean = Classimg_data.updateMask(Orchard_ClassifiedImg).sample(
  region = region,
  scale = 10,
  numPixels = 1000,
  seed = 42,
  tileScale = 16)

