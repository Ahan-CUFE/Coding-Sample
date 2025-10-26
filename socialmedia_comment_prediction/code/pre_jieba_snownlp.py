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
        saved=['a','v','n','vn']




        cutwordscopy = [x for x in cutwords if x[1] in saved]
        cutwordscopy=list(map(str,cutwords))
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
datasetnew.to_csv('/Users/gaoyihan/pythonProject/socialmedia_comment_prediction/pre_jieba_snownlp_avn_vn_saved.csv')


##########################################
#目前来看，snownlp的打分效果不是十分理想，但是足以进行初步改进，后续模型改进有两个思路
#（1）增添积极/消极 语料库 利用sonwnlp的训练模块进行模型的适应性训练。
#（2）使用finbert包再进行建模。
#