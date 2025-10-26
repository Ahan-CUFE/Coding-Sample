import pandas as pd
from snownlp import SnowNLP

#Only need to run once to train the model
from snownlp import sentiment
sentiment.train('/Users/gaoyihan/opt/anaconda3/envs/nlp/lib/python3.8/site-packages/snownlp/sentiment/negtive.txt', '/Users/gaoyihan/opt/anaconda3/envs/nlp/lib/python3.8/site-packages/snownlp/sentiment/positive.txt')
sentiment.save('fintext.marshal')


dataset = pd.read_excel('/Users/gaoyihan/Desktop/新沪深300股评文本.xlsx')
from opencc import OpenCC

#transfer the traditional Chinese to simplified Chinese

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
#non-text information filtering
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
# # keep adjectives
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

        ####Keep some necessary parts of speech for text analysis.
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




