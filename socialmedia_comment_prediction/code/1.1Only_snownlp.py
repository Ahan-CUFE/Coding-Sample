import pandas as pd
import jieba
from snownlp import SnowNLP
dataset = pd.read_csv('/Users/gaoyihan/Desktop/hs300text.csv')

dataset['emotionvalue'] = 0

for i in range(len(dataset.emotionvalue)):
    ss=SnowNLP(dataset.iloc[i,0]).sentiments
    dataset.loc[i,'emotionvalue'] = ss

dataset.to_csv('/Users/gaoyihan/pythonProject/socialmedia_comment_prediction/onlySnownlp.csv')