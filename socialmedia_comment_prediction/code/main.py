import pandas as pd
import matplotlib.pyplot as plt
from pandas import Series

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

onlysnowline=onlysnow.resample('1D',on='date')['emotionvalue'].mean()
preandsnow =preandsnow.resample('1D',on = 'date')['emotionvalue'].mean()
prejiebasnow = prejiebasnow.resample('1D',on = 'date')['emotionvalue'].mean()


onlysnowline.plot()
preandsnow.plot()
prejiebasnow.plot()
plt.show()

#####################################################################
#导入hs300数据

