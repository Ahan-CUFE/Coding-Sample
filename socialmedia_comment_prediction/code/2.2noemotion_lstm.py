import pandas as pd
import matplotlib.pyplot as plt
from pandas import Series
from factor_analyzer.factor_analyzer import calculate_bartlett_sphericity

dataset = pd.read_csv('/Users/gaoyihan/pythonProject/socialmedia_comment_prediction/advantaged_trained_nlp.csv',index_col=0)
dataset['date'] =pd.to_datetime(dataset.发帖时间)
dataset = dataset.resample('1D',on = 'date')['emotionvalue'].mean()


datasetnew = {'date':dataset.index,'emotionvalue':dataset.values}
datasetnew = pd.DataFrame(datasetnew)
datasetnew['date'] = pd.to_datetime(datasetnew.date).dt.date
datasetnew.set_index('date',inplace=True)

#Import hs300 data to train fundamental factor prediction model.
hs300 = pd.read_csv('/Users/gaoyihan/Desktop/hs300.csv')
hs300['date'] = pd.to_datetime(hs300.date).dt.date
# hs300['date1']=hs300.date
hs300.set_index('date',inplace=True)
# hs300['emotionvalue'] = 0
#
# for i in hs300.index:
#     hs300.loc[i,'emotionvalue'] = datasetnew.loc[i,'emotionvalue']

# hs300=hs300.drop('成交额')

from sklearn import preprocessing
hs300.fillna(0,inplace=True)

data_scaled = (hs300 - hs300.min()) / (hs300.max()-hs300.min())



# PCA 

chi_square_value, p_value = calculate_bartlett_sphericity(hs300)
print(chi_square_value, p_value)
# KMO

from factor_analyzer.factor_analyzer import calculate_kmo

kmo_all, kmo_model = calculate_kmo(hs300)
print(kmo_all)

print(kmo_all)
from sklearn.decomposition import PCA
# pca = PCA(n_components=4)
pca = PCA(n_components=6)
pca.fit(data_scaled)
pca.explained_variance_ratio_
datapca= pca.transform(data_scaled)
datapca = pd.DataFrame(datapca)
# 0.60298879, 0.13669649, 0.09426786, 0.06602151, 0.03756129,0.03062302] = 0.94

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


traindata = datapca[:-int((0.2*len(datapca.index)))]
testdata  = datapca[-int((0.2*len(datapca.index))):]

traindataY=data_scaled[:-int((0.2*len(datapca.index)))]
testdataY  = data_scaled[-int((0.2*len(datapca.index))):]

print('traindata shape is ',traindata.shape,'testdata shape is ',testdata.shape)



def creatXY(datasetX,datasetY,n_past):
    dataX=[]
    dataY=[]
    for i in range(n_past,len(datasetX)):
        dataX.append(datasetX.iloc[i-n_past:i,0:datasetX.shape[1]])
        dataY.append(datasetY.iloc[i,0])  # 把次日开盘价作为预测标签

    return np.array(dataX),np.array(dataY)
# traindataYnp = traindataY.to_numpy()
# testdataYnp  = testdataY.to_numpy()

# traindataY.reset_index(inplace=True,drop=True)
# testdataY.reset_index(inplace=True,drop=True)

n_past=7
trainX,trainY=creatXY(traindata,traindataY,n_past)
testX,testY = creatXY(testdata,testdataY,n_past)


print("trainX Shape-- ",trainX.shape)
print("trainY Shape-- ",trainY.shape)

print("testX Shape-- ",testX.shape)
print("testY Shape-- ",testY.shape)


def model_creator(optimizer):
    grid_model=Sequential()
    grid_model.add(LSTM(50,return_sequences=True,input_shape=(n_past,6)))
    grid_model.add(LSTM(50))
    grid_model.add(Dropout(0.2))
    grid_model.add(Dense(1))

    grid_model.compile(loss = 'mse',optimizer=optimizer)
    return grid_model
grid_model=KerasRegressor(build_fn=model_creator('adam'),verbose=1)
grid_model.fit(trainX, trainY, epochs=300, batch_size=64, verbose=1)
Y_hat = grid_model.predict(testX)
Y_hat_unscaled = (hs300.open.max()-hs300.open.min())*Y_hat + hs300.open.min()
Y_test_unscaled = (hs300.open.max()-hs300.open.min())*testY + hs300.open.min()


plt.plot(Y_hat_unscaled,label ='y_hat')
plt.plot(Y_test_unscaled,label='t_true')
plt.legend()
plt.savefig('/Users/gaoyihan/Desktop/结果/only_past7_pre1_epoch300.jpg')
plt.show()

def calculate_mse(y_hat,y_true):
    detla = y_hat -y_true
    mse = 0
    for i in range(len(detla)):
        mse = mse+detla[i] *detla[i]

    mse = mse / len(detla)
    return mse
mse = calculate_mse(Y_hat_unscaled,Y_test_unscaled)
print(mse)


