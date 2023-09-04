
import geemap
import ee
geemap.set_proxy(port=)
ee.Initialize()
Map = geemap.Map()
Map
    
# 输入参数
Roi = 'users/wanghan02191127/shanxilinyi'
ImageYear = 2021
savename = 'Shaanxilinyi_Other1AndOrhchard2_ClassifiedImg'
samp =' users/wanghan02191127/shanxilinyi-sample'

dapeng = samp.filter(ee.Filter.eq('Type', 1)) \
            .map(lambda fea: fea.set('class', 0))
li = samp.filter(ee.Filter.eq('Type', 2)) \
         .map(lambda fea: fea.set('class', 1))
pingguo = samp.filter(ee.Filter.eq('Type', 3)) \
             .map(lambda fea: fea.set('class', 2))
shiliu = samp.filter(ee.Filter.eq('Type', 4)) \
            .map(lambda fea: fea.set('class', 3))
shizi = samp.filter(ee.Filter.eq('Type', 5)) \
           .map(lambda fea: fea.set('class', 4))
tao = samp.filter(ee.Filter.eq('Type', 6)) \
         .map(lambda fea: fea.set('class', 5))
yumi = samp.filter(ee.Filter.eq('Type', 7)) \
          .map(lambda fea: fea.set('class', 6))
zao = samp.filter(ee.Filter.eq('Type', 8)) \
         .map(lambda fea: fea.set('class', 7))

# applepoint_style = {'color': 'red'}
# otherpoint_style = {'color': 'green'}
# Map.addLayer(pingguo,applepoint_style,'applepoint')
# Map.addLayer(dapeng,otherpoint_style,'otherpoint')
# 样本点与研究区
# 样本点OrchardAndOtherClass：果园为1，非果园为0
# print("appletNumber",pingguo.size())
# print("dapengNumber",dapeng.size())

# 研究区边界显示
# Roi_outline=ee.Image().toByte().paint(featureCollection = Roi,color = 0,width = 0.5)
# roi_style = {'palette': 'red'}
# Map.addLayer(Roi_outline,roi_style,'Roi_outline')
# Map.centerObject(Roi,6)

# 已有土地利用数据
imageVisParam_ESA = {"opacity":1,
                     "bands":["Map"],
                     "min":10,
                     "max":110,
                    "palette":["006400","ffbb22","ffff4c","f096ff","fa0000","b4b4b4","f0f0f0","0064c8","0096a0","00cf75","fae6a0"]
                    }
ESA_Roi = ee.ImageCollection("ESA/WorldCover/v100")\
            .filterBounds(Roi)\
            .select('Map')\
            .mosaic()\
            .clip(Roi)
crop_Roi = ESA_Roi.eq(ee.Image.constant(40))  
crop_Roi_mask = crop_Roi.updateMask(crop_Roi.mask())
Map.addLayer(ESA_Roi,imageVisParam_ESA,'ESA_WorldCover')
Map.addLayer(crop_Roi_mask, {'min': 0, 'max': 1, 'palette': ['black', 'green']}, 'ESA_Cropland')

# 哨兵影像数据
year = ee.Number(ImageYear)
startDay = ee.Date.fromYMD(year,1,1)
endDay = ee.Date.fromYMD(year,12,30)
bands =ee.list(['B2','B3','B4','B5','B6','B7','B8','B11','B12','NDVI','NDRE1','A','EVI','RESI','NDBI'])

# Function
# S2 ReProcess Function
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
    ndre1 = image.normalizedDifference(['B6', 'B5'])
    resi = image.expression('(RE3+RE2-RE1)/(RE3+RE2+RE1)', {
        'RE1': image.select('B5'),
        'RE2': image.select('B6'),
        'RE3': image.select('B7')
    })
    A = image.expression('(RE3-RE2)/RE1', {
        'RE1': image.select('B5'),
        'RE2': image.select('B6'),
        'RE3': image.select('B7')
    })
    return image.addBands(ndvi.rename("NDVI")) \
                .addBands(ndre1.rename('NDRE1')) \
                .addBands(resi.rename("RESI")) \
                .addBands(A.rename("A")) 

#  全年有效影像——原始波段以及计算波段
Sentinel2 = ee.ImageCollection("COPERNICUS/S2_SR")\
              .filterDate(startDay,endDay)\
              .filterBounds(Roi)\
              .map(sentinel2toa)\
              .map(cloudMask)\
              .map(addIndices)
S2_ValidiImg_clipRoi = Sentinel2.select(bands)\
                                .map(lambda image:image.clip(Roi))
                                
Map.addLayer(S2_ValidiImg_clipRoi.median(),{min: 0, max: 0.3, bands: ['B4', 'B3', 'B2']},'S2_ValidiImg_clipRoi_RGB')

# 2.月均值合成影像——Bands
month = ee.List.sequence(1,12,1)
Month_Com = []
for i in month:
    i = ee.Number(i)
    SDay = ee.Date.fromYMD(year, i, 1)
    EDay = SDay.advance(1, 'month')
    Img = S2_ValidiImg_clipRoi.filterDate(SDay, EDay).mean()
    Img = Img.set('time', SDay.format('yyyy-MM-dd'))
    Month_Com.append(Img)
    
    

# 3.NDVI曲线谐波回归平滑
# harmonic regression 谐波回归
# c0 + c1*cos(2*pi*t) + c2*sin(2*pi*t) = NDVI
time_reference = startDay
def process_image(img):
    tstamp_NDVI = ee.Date(img.get('time'))
    tdelta_NDVI = tstamp_NDVI.difference(time_reference, 'year')
    # 构建用于拟合方程的图像
    img_fitting_NDVI = img.select() \
        .addBands(1)\
        .addBands(ee.Image(tdelta_NDVI.multiply(2 * 3.141592653589793).cos()).rename('cos')) \
        .addBands(ee.Image(tdelta_NDVI.multiply(2 * 3.141592653589793).sin()).rename('sin')) \
        .addBands(img.select('NDVI')) \
        .toDouble()
    return img_fitting_NDVI
s2_Month_Com_fitting_NDVI = Month_Com.map(process_image)

dependent_NDVI = ee.List(['NDVI'])
harmonicIndependents_NDVI= ee.List(['constant', 'cos', 'sin'])
# The output of the regerssion reduction is a[X,Y] array image.
harmonic_NDVI = s2_Month_Com_fitting_NDVI.select(harmonicIndependents_NDVI.cat(dependent_NDVI)) \
                                          .reduce(ee.Reducer.linearRegression(
                                              harmonicIndependents_NDVI.length(), dependent_NDVI.length()))
coefficients_NDVI = harmonic_NDVI.select('coefficients')\
                                 .matrixTranspose() \
                                 .arrayFlatten([dependent_NDVI, harmonicIndependents_NDVI])\
                                 .clip(Roi)  

def add_fitted_NDVI(image):
    fitted_NDVI = image.select(harmonicIndependents_NDVI) \
                        .multiply(coefficients_NDVI) \
                        .reduce('sum') \
                        .rename('fitted_NDVI')
    return image.addBands(fitted_NDVI)
s2_Month_Com_fitted_NDVI = s2_Month_Com_fitting_NDVI.map(add_fitted_NDVI)

Map.addLayer(s2_Month_Com_fitted_NDVI.select(['NDVI','fitted_NDVI']),{},'s2_Month_Com_fitted_NDVI')

# 4.NDRE1曲线谐波回归平滑

time_reference = startDay
def process_image(img):
    tstamp_NDRE1 = ee.Date(img.get('time'))
    tdelta_NDRE1 = tstamp_NDRE1.difference(time_reference, 'year')
    # 构建用于拟合方程的图像
    img_fitting_NDRE1 = img.select() \
        .addBands(1)\
        .addBands(ee.Image(tdelta_NDRE1.multiply(2 * 3.141592653589793).cos()).rename('cos')) \
        .addBands(ee.Image(tdelta_NDRE1.multiply(2 * 3.141592653589793).sin()).rename('sin')) \
        .addBands(img.select('NDRE1')) \
        .toDouble()
    return img_fitting_NDRE1
s2_Month_Com_fitting_NDRE1 = Month_Com.map(process_image)

dependent_NDRE1 = ee.List(['NDRE1'])
harmonicIndependents_NDRE1= ee.List(['constant', 'cos', 'sin'])
# The output of the regerssion reduction is a[X,Y] array image.
harmonic_NDRE1 = s2_Month_Com_fitting_NDRE1.select(harmonicIndependents_NDRE1.cat(dependent_NDRE1)) \
                                          .reduce(ee.Reducer.linearRegression(
                                              harmonicIndependents_NDRE1.length(), dependent_NDRE1.length()))
coefficients_NDRE1 = harmonic_NDRE1.select('coefficients')\
                                 .matrixTranspose() \
                                 .arrayFlatten([dependent_NDRE1, harmonicIndependents_NDRE1])\
                                 .clip(Roi)  

def add_fitted_NDRE1(image):
    fitted_NDRE1 = image.select(harmonicIndependents_NDRE1) \
                        .multiply(coefficients_NDRE1) \
                        .reduce('sum') \
                        .rename('fitted_NDRE1')
    return image.addBands(fitted_NDRE1)
s2_Month_Com_fitted_NDRE1 = s2_Month_Com_fitting_NDRE1.map(add_fitted_NDRE1)

Map.addLayer(s2_Month_Com_fitted_NDRE1.select(['NDRE1','fitted_NDRE1']),{},'s2_Month_Com_fitted_NDRE1')


# 5.A曲线谐波回归平滑
time_reference = startDay
def process_image(img):
    tstamp_A = ee.Date(img.get('time'))
    tdelta_A = tstamp_A.difference(time_reference, 'year')
    # 构建用于拟合方程的图像
    img_fitting_A = img.select() \
        .addBands(1)\
        .addBands(ee.Image(tdelta_A.multiply(2 * 3.141592653589793).cos()).rename('cos')) \
        .addBands(ee.Image(tdelta_A.multiply(2 * 3.141592653589793).sin()).rename('sin')) \
        .addBands(img.select('A')) \
        .toDouble()
    return img_fitting_A
s2_Month_Com_fitting_A = Month_Com.map(process_image)

dependent_A = ee.List(['A'])
harmonicIndependents_A= ee.List(['constant', 'cos', 'sin'])
# The output of the regerssion reduction is a[X,Y] array image.
harmonic_A = s2_Month_Com_fitting_A.select(harmonicIndependents_A.cat(dependent_A)) \
                                          .reduce(ee.Reducer.linearRegression(
                                              harmonicIndependents_A.length(), dependent_A.length()))
coefficients_A = harmonic_A.select('coefficients')\
                                 .matrixTranspose() \
                                 .arrayFlatten([dependent_A, harmonicIndependents_A])\
                                 .clip(Roi)  

def add_fitted_A(image):
    fitted_A = image.select(harmonicIndependents_A) \
                        .multiply(coefficients_A) \
                        .reduce('sum') \
                        .rename('fitted_A')
    return image.addBands(fitted_A)
s2_Month_Com_fitted_A = s2_Month_Com_fitting_A.map(add_fitted_A)

Map.addLayer(s2_Month_Com_fitted_A.select(['A','fitted_A']),{},'s2_Month_Com_fitted_A')

# 6.纹理特征
S2_5Month = S2_ValidiImg_clipRoi.filterDate('2021-06-01','2021-11-01').median()
S2_5Month_ndvi =S2_5Month.clip(Roi).select('NDVI')
S2_5Month_ndvi_Uint16 = S2_5Month_ndvi.add(127.5)\
                                      .multiply(127.5)\
                                      .toUint16()
S2_5Month_ndvi_glcm = S2_5Month_ndvi_Uint16.glcmTexture()


#7.DEM特征
dataset = ee.image('USGS/SRTMGL1_003')
elevation = dataset.select('elevation')
slope = ee.Terrain.slope(elevation)
hillshade = ee.Terrain.hillshade(elevation)


# 8.分类
def mergebands(current, previous):
    return ee.Image(previous).addBands(current)

S2_Month_Com_NDVI = ee.Image(s2_Month_Com_fitted_NDVI.select(['fitted_NDVI'])
                             .iterate(mergebands, ee.Image([])))
S2_Month_Com_NDRE1 = ee.Image(s2_Month_Com_fitted_NDRE1.select(['fitted_NDRE1'])
                              .iterate(mergebands, ee.Image([])))
S2_Month_Com_A = ee.Image(s2_Month_Com_fitted_A.select(['fitted_A'])
                          .iterate(mergebands, ee.Image([])))

S2_Median_Year = S2_ValidiImg_clipRoi.median()

glcmb=['NDVI_asm','NDVI_contrast','NDVI_corr','NDVI_ent','NDVI_idm']

Pre_snic_Image = ee.Image.cat([S2_Month_Com_NDVI,
                                 S2_Month_Com_NDRE1,
                                 S2_Month_Com_A,
                                 S2_Median_Year,
                                 S2_5Month_ndvi_glcm.select(glcmb),
                                 elevation,
                                 slope,
                                 hillshade
                                 ]).updateMask(crop_Roi)


print(Pre_snic_Image)




seeds = ee.Algorithms.Image.Segmentation.seedGrid(25)
ClassInPut_Image = ee.Algorithms.Image.Segmentation.SNIC(image=Pre_snic_Image,
                                            size=32,
                                            compactness=0,
                                            connectivity=8,
                                            neighborhoodSize=256,
                                            seeds=seeds)
# 先分地块样本在提取像素值，验证
print('先分地点，在提取像素值，验证')

SampleFea_rondom = samp.randomColumn('random',42)
tariningPoints = SampleFea_rondom.filter(ee.Filter.lt('random',0.7))
testingPoints = SampleFea_rondom.filter(ee.Filter.gte('random',0.7))
trainingPoints_value = ClassInPut_Image.sampleRegions(
    collection = tariningPoints,
    properties = ['Type'],
    scale = 10,
    tileScale = 16,
    geometries = True
)
testingPoints_value = ClassInPut_Image.sampleRegions(
    collection = testingPoints,
    properties = ['Type'],
    scale = 10,
    tileScale =16,
    geometries = True
)
# print(trainingPoints_value.size(),testingPoints_value.size())

# 分类
InputImgName = ClassInPut_Image.bandNames()

TrainedModel1 = ee.Classifier.smileRandomForest(
    numberOfTrees = 100,
    seed = 42
).tarin(trainingPoints_value,'Type',InputImgName)

ClassifiedResult1 = ClassInPut_Image.select(InputImgName).classify(TrainedModel1)
classified_result_vis_params = {
    'min': 0,
    'max': 1,
    'palette': ['green', 'red']
}
Map.addLayer(ClassifiedResult1, classified_result_vis_params, 'Classified Result')

ValidatedResult1 = testingPoints_value.classify(TrainedModel1)
ValidatedAccuracy1 = ValidatedResult1.errorMatrix('Type', 'classification')

print('ErrorMatrix:',ValidatedAccuracy1)
print('Validate accuracy:',ValidatedAccuracy1.accuracy())
print('User acc:',ValidatedAccuracy1.consumersAccuracy())
print('Prod acc:',ValidatedAccuracy1.producersAccuracy())
print('Kappa:',ValidatedAccuracy1.kappa())

# 导出影像
ClassifiedImg = ClassifiedResult1.toByte()
ClassifiedImg = ClassifiedImg.remap([0,1], [1,2])

# Define export parameters
export_image_args = {
    'image': ClassifiedImg,
    'description': savename,
    'region': Roi,
    'scale': 10,
    'maxPixels': 1e13
}

# Export to Google Drive
geemap.ee_export_image_to_drive(**export_image_args, folder='LandUsePaper')

# Export to Earth Engine Asset
geemap.ee_export_image_to_asset(**export_image_args)










