import pandas as pd
import matplotlib.pyplot as plt
from pandas import Series
from factor_analyzer.factor_analyzer import calculate_bartlett_sphericity
onlysnow = pd.read_csv('/Users/gaoyihan/pythonProject/socialmedia_comment_prediction/onlySnownlp.csv',index_col=0)
preandsnow = pd.read_csv('/Users/gaoyihan/pythonProject/socialmedia_comment_prediction/processed_text_and_snownlp.csv',index_col=0)
prejiebasnow = pd.read_csv('/Users/gaoyihan/pythonProject/socialmedia_comment_prediction/pre_jieba_snownlp_avn_vn_saved.csv',index_col=0)
value1 = onlysnow.emotionvalue
value2 = preandsnow.emotionvalue
value3 = prejiebasnow.emotionvalue

onlysnow.head()

onlysnow['date'] =pd.to_datetime(onlysnow.发帖时间)
preandsnow['date'] = pd.to_datetime(preandsnow.发帖时间)
prejiebasnow['date']=pd.to_datetime(prejiebasnow.发帖时间)

# onlysnowline=onlysnow.resample('1D',on='date')['emotionvalue'].mean()
preandsnow =preandsnow.resample('1D',on = 'date')['emotionvalue'].mean()
prejiebasnow = prejiebasnow.resample('1D',on = 'date')['emotionvalue'].mean()


# onlysnowline.plot()
# preandsnow.plot()
# prejiebasnow.plot()
# plt.show()

#####################################################################
#导入hs300数据,训练基本面因子的预测模型。
hs300 = pd.read_csv('/Users/gaoyihan/pythonProject/socialmedia_comment_prediction/hs300data.csv')
hs300['date'] = pd.to_datetime(hs300.date)
hs300.set_index('date',inplace=True)


#PCA 降维  ，
# 球状检验，即检验各个变量是否各自独立。
chi_square_value, p_value = calculate_bartlett_sphericity(hs300)
print(chi_square_value, p_value)
# KMO检验
# 检查变量间的相关性和偏相关性，取值在0-1之间；KOM统计量越接近1，变量间的相关性越强，偏相关性越弱，因子分析的效果越好。
# 通常取值从0.6开始进行因子分析
from factor_analyzer.factor_analyzer import calculate_kmo

kmo_all, kmo_model = calculate_kmo(hs300)
print(kmo_all)
from sklearn.decomposition import PCA
# pca = PCA(n_components=4)
pca = PCA(n_components=2)

## when n_components = 4 ,explained_variance_ratio_=[0.80672961, 0.14377399, 0.03827036, 0.01122587]
# so just choose n=2 can save the most of the information of the factor
pca.fit(hs300)
# hs300= pca.transform(hs300)
# hs300 = pd.DataFrame(hs300)
# # pca.explained_variance_ratio_
# hs300_d = pd.read_csv('/Users/gaoyihan/pythonProject/socialmedia_comment_prediction/hs300data.csv')
#
# hs300['date'] = hs300_d['date']
# hs300['date'] = pd.to_datetime(hs300.date)
# hs300.set_index('date',inplace=True)




import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from keras.models import Sequential
from keras.layers import LSTM
from keras.layers import Dense, Dropout
from keras import optimizers
from sklearn.model_selection import GridSearchCV
from scikeras.wrappers import KerasRegressor
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'


traindata = hs300[:-28]
testdata  = hs300[-28:]

print('traindata shape is ',traindata.shape,'testdata shape is ',testdata.shape)


##数据归一化，增快优化算法收敛速度

# 将数据scale为(0,1)间数据
traindata_scaled = (traindata - traindata.min()) / (traindata.max()-traindata.min())
testdata_scaled = (testdata - testdata.min()) / (testdata.max()-testdata.min())


### 创建 属性 和 标签 的函数:dataset为输入的数据集，npast为滞后的步数。
def creatXY(dataset,n_past):
    dataX=[]
    dataY=[]
    for i in range(n_past+1,len(dataset)):
        dataX.append(dataset.iloc[i-n_past-1:i-1,0:dataset.shape[1]])
        dataY.append(dataset.iloc[i,0])  # 把次日开盘价作为预测标签

    return np.array(dataX),np.array(dataY)

trainX,trainY=creatXY(traindata_scaled,5)
testX,testY = creatXY(testdata_scaled,5)

print("trainX Shape-- ",trainX.shape)
print("trainY Shape-- ",trainY.shape)

print("testX Shape-- ",testX.shape)
print("testY Shape-- ",testY.shape)



def model_creator(optimizer):
    grid_model=Sequential()
    grid_model.add(LSTM(50,return_sequences=True,input_shape=(5,13)))
    grid_model.add(LSTM(50))
    grid_model.add(Dropout(0.2))
    grid_model.add(Dense(1))

    grid_model.compile(loss = 'mse',optimizer=optimizer)
    return grid_model

grid_model=KerasRegressor(build_fn=model_creator('adam'),verbose=1)
##超参数设置
# parameters = {'batch_size' : [16,20],
#             'epochs' : [8,10],
# }
# grid_search = GridSearchCV(estimator = grid_model,param_grid = parameters , cv = 2)
# grid_search = grid_search.fit(trainX,trainY)

grid_model.fit(trainX, trainY, epochs=100, batch_size=64, verbose=1)

Y_hat = grid_model.predict(testX)

Y_hat_unscaled = (testdata.open.max()-testdata.open.min())*Y_hat + testdata.open.min()
Y_test_unscaled = (testdata.open.max()-testdata.open.min())*testY + testdata.open.min()


plt.plot(Y_hat_unscaled)
plt.plot(Y_test_unscaled)
plt.show()



#####最佳超参数网格化搜索，优化器等等，网络层数。
#####样本数增加至500条
#####增加数据维度，增加解释变量的种类：速动比率，流通速率等等


