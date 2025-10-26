import pandas as pd
from opencc import OpenCC
from snownlp import SnowNLP

dataset = pd.read_csv('/Users/gaoyihan/Desktop/hs300text.csv')
from opencc import OpenCC

#繁体字转换为简体字

def chinese_standard(text:str, conversion='t2s'):
    cc = OpenCC(conversion)
    return cc.convert(text)

#make a copy of dataset to store the new dataset with the processed text

datasetnew= dataset.copy()
for i in range(len(dataset['标题'])):
    datasetnew.loc[i,'标题'] = chinese_standard(dataset.loc[i,'标题'])

#非文本信息过滤
import re
def clear_character(text):
    pattern = [
        "[^\u4e00-\u9fa5^a-z^A-Z^0-9^\u0020^\u0027^\u002e]",  # save_standing_character
        "\.$"
    ]
    return re.sub('|'.join(pattern),'',text)

for i in range(len(dataset['标题'])):
    datasetnew.loc[i,'标题'] = clear_character(datasetnew.loc[i,'标题'])


#calcute the sccore of the new- processed text
#some

datasetnew['emotionvalue'] = 0
datasetnew=datasetnew.dropna()
for i in range(len(datasetnew.emotionvalue)):
    try:
        ss=SnowNLP(datasetnew.iloc[i,0]).sentiments
        datasetnew.loc[i,'emotionvalue'] = ss
    except :
        print('there is something wrong with index',i)


datasetnew.to_csv('/Users/gaoyihan/pythonProject/socialmedia_comment_prediction/processed_text_and_snownlp.csv')
