import pandas as pd
from snownlp import SnowNLP

#只需要train一次
# from snownlp import sentiment
# sentiment.train('/Users/gaoyihan/opt/anaconda3/envs/nlp/lib/python3.8/site-packages/snownlp/sentiment/negtive.txt', '/Users/gaoyihan/opt/anaconda3/envs/nlp/lib/python3.8/site-packages/snownlp/sentiment/positive.txt')
# sentiment.save('fintext.marshal')


dataset = pd.read_excel('/Users/gaoyihan/Desktop/新沪深300股评文本.xlsx')
from opencc import OpenCC

#繁体字转换为简体字

def chinese_standard(text:str, conversion='t2s'):
    cc = OpenCC(conversion)
    return cc.convert(text)

#make a copy of dataset to store the new dataset with the processed text

datasetnew= dataset.copy()
for i in range(len(datasetnew['标题'])):
    try:
        datasetnew.loc[i,'标题'] = chinese_standard(datasetnew.loc[i,'标题'])
    except:
        print('something wrong with ',i)
        datasetnew.drop(i,axis=0,inplace=True)
datasetnew = datasetnew.reset_index(drop=True)
#非文本信息过滤
import re
def clear_character(text):
    pattern = [
        "[^\u4e00-\u9fa5^a-z^A-Z^0-9^\u0020^\u0027^\u002e]",  # save_standing_character
        "\.$"
    ]
    return re.sub('|'.join(pattern),'',text)

for i in range(len(datasetnew['标题'])):
    try:
        datasetnew.loc[i,'标题'] = clear_character(datasetnew.loc[i,'标题'])
    except:
        print('something wrong with ',i)

#calcute the sccore of the new- processed text
#considering that the Snownlp don't doing well in Chinese Word Segmentation,so we will use jieba to do
#Chinese word Seqmentation ,and sooner use Snownlp to scorre that
import jieba
import jieba.posseg as psg

#
# words5 = [ (w.word, w.flag) for w in psg.cut(words1) ]
# # 保留形容词
# saved = ['a',]
# words5 =[x for x in words5 if x[1] in saved]
# print(words5)
# # [('快', 'a'), ('好', 'a'), ('好', 'a'), ('不错', 'a')]

datasetnew['emotionvalue'] = 0
datasetnew['cuttedwordslist'] = '无'
for i in range(len(datasetnew.emotionvalue)):
    try:
        cuttedobj = datasetnew.loc[i,'标题']
        cutwords = [(w.word,w.flag) for w in psg.cut(cuttedobj)]

        ####酌情考虑保存一些必要词性做文本分析。
        saved=['a','v','n','vn','d','vd','i']

        # cutwordscopy=list(map(str,cutwords))
        cutwordscopy = [x[0] for x in cutwords if x[1] in saved]
        cutwordscopy = ''.join(cutwordscopy)
        datasetnew.loc[i,'cuttedwordslist']= cutwordscopy  #to recoard


###calcute the sccore
        ss = SnowNLP(cutwordscopy)
        datasetnew.loc[i,'emotionvalue'] = ss.sentiments
    except:
        print('there is something wrong with index',i)

#
# for i in range(len(datasetnew.emotionvalue)):
#     try:
#         ss=SnowNLP(datasetnew.iloc[i,0]).sentiments
#         datasetnew.loc[i,'emotionvalue'] = ss
#     except :
#         print('there is something wrong with index',i)
datasetnew.to_csv('/Users/gaoyihan/pythonProject/socialmedia_comment_prediction/nnadvantaged_trained_nlp.csv')


###### train nlp modle with fin_text






##########################################
#目前来看，snownlp的打分效果不是十分理想，但是足以进行初步改进，后续模型改进有两个思路
#（1）增添积极/消极 语料库 利用sonwnlp的训练模块进行模型的适应性训练。
#（2）使用finbert包再进行建模。
#