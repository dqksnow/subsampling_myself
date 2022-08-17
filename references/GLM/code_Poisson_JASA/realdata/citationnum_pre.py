# -*- coding: utf-8 -*-
"""
Created on Sat Jan 12 14:23:51 2019

@author: yujun
"""
import json
import pandas as pd
import os
#from fuzzywuzzy import fuzz
#from fuzzywuzzy import process
import numpy as np

f_j = open(r'D:/讨论班/report5/data/sjrclean_2017.csv',encoding='utf8')
journal_info = pd.read_csv(f_j,sep=';')


def journalcheck(venue):
    n=journal_info.shape[0]
    for i in range(n):
        title = journal_info.loc[i]['Title'].lower()
        title = title.replace("&","and")
        pattern = venue.lower()
        pattern = pattern.replace("&","and")
        if title == pattern:
            journal = journal_info.loc[i]
            res = []
            n = len(journal)
            for j in range(1, n):
                try:
                    res.append(journal[j])
                except:
                    res.append('NA')
            return res
    for i in range(n):
        title = journal_info.loc[i]['Title'].lower()
        title = title.replace("&","and")
        if title.find(pattern) != -1:
             journal = journal_info.loc[i]
             res = []
             n = len(journal)            
             for j in range(1, n):
                try:
                     res.append(journal[j])
                except:
                     res.append('NA')
             return res
    return ['NA','NA','NA','NA']

journaltitle=journal_info.loc[:]['Title'].tolist()
journaltitle1 = [journaltitle[i].lower().replace('&',"and") for i in range(len(journaltitle))]
def journalchecknew(venue,year):
    temp1 = [i for i,x in enumerate(journaltitle1) if x.find(venue)!=-1]   
    if temp1==[]:
        return ['NA']
    temp3 = [i for i,x in enumerate(journaltitle1) if x==venue]
    if temp3!=[]:
        return journaltitle1[temp3[0]]
    temp2 = [i for i,x in enumerate(journaltitle1) if x.find(venue)!=-1&x.find(str(year))!=-1]
    if temp2!=[]:
        return journaltitle1[temp2[0]]
    else:
        return journaltitle1[temp1[0]]

b = pd.DataFrame(journaltitle1,columns=['journaltitle1'])
sjrmerge= pd.concat([journal_info,b],axis=1)
sjrmerge.to_csv('D:/讨论班/report5/data/sjrmerge_2017.csv')



datadir = r'D:/讨论班/report5/data/dblp.v10/dblp-ref/dblp-ref-0.json'
f=open(datadir,'r')
a=f.readlines()        
authordict= dict()
authorimpactdict= dict()
for i in range(len(a)):
    aritcle = json.loads(a[i])
    for author in aritcle['authors']:
        if authordict.get(author):
            authordict[author] += 1
            authorimpactdict[author] += aritcle['n_citation']
        else:
            authordict[author] = 1
            authorimpactdict[author] = aritcle['n_citation']
print("finish"+datadir)

datadir = r'D:/讨论班/report5/data/dblp.v10/dblp-ref/dblp-ref-1.json'
f=open(datadir,'r')
a=f.readlines()        
for i in range(len(a)):
    aritcle = json.loads(a[i])
    if 'authors' in aritcle.keys():
        for author in aritcle['authors']:
            if authordict.get(author):
                authordict[author] += 1
                authorimpactdict[author] += aritcle['n_citation']
            else:
                authordict[author] = 1
                authorimpactdict[author] = aritcle['n_citation']
    else:
        continue
print("finish"+datadir)
        
datadir = r'D:/讨论班/report5/data/dblp.v10/dblp-ref/dblp-ref-2.json'
f=open(datadir,'r')
a=f.readlines()        
for i in range(len(a)):
    aritcle = json.loads(a[i])
    if 'authors' in aritcle.keys():
        for author in aritcle['authors']:
            if authordict.get(author):
                authordict[author] += 1
                authorimpactdict[author] += aritcle['n_citation']
            else:
                authordict[author] = 1
                authorimpactdict[author] = aritcle['n_citation']
    else:
        continue
print("finish"+datadir)
        
                      
datadir = r'D:/讨论班/report5/data/dblp.v10/dblp-ref/dblp-ref-3.json'
f=open(datadir,'r')
a=f.readlines()        
for i in range(len(a)):
    aritcle = json.loads(a[i])
    if 'authors' in aritcle.keys():
        for author in aritcle['authors']:
            if authordict.get(author):
                authordict[author] += 1
                authorimpactdict[author] += aritcle['n_citation']
            else:
                authordict[author] = 1
                authorimpactdict[author] = aritcle['n_citation']
    else:
        continue
print("finish"+datadir)        
                    


datafinal = pd.DataFrame()
datadir = r'D:/讨论班/report5/data/dblp.v10/dblp-ref/dblp-ref-0.json'
f=open(datadir,'r')
a=f.readlines()
for i in range(len(a)):
    aritcle = json.loads(a[i])
    if 'abstract' in aritcle.keys():
        abs_words= len(aritcle['abstract'].split())
    else:
        abs_words= 0
    if 'authors' in aritcle.keys():
        number_authors = len(aritcle['authors'])
    else:
        continue
    citation = aritcle['n_citation']
    if 'references' in aritcle.keys():
        number_ref = len(aritcle['references'])
    else:
        number_ref= 0 
    title = aritcle['title']
    num_title = len(aritcle['title'].split())
    venue = aritcle['venue'].lower().replace("&","and")
    year = aritcle['year']
    venue1 = journalchecknew(venue,year)
    authorimacttemp = []
    authorpapertemp = []
    for author in aritcle['authors']:
        authorimacttemp.append(authorimpactdict[author])
        authorpapertemp.append(authordict[author])
    authorimapct1=max(authorimacttemp)
    authorimapct2=np.mean(authorimacttemp)
    authorpaper1=max(authorpapertemp)
    authorpaper2=np.mean(authorpapertemp)
        
    #venuematch = process.extractOne(venue,journaltitle)
#    journal_info0 = journalcheck(venue)[0]
#    journal_info1 = journalcheck(venue)[1]
#    journal_info2 = journalcheck(venue)[2]
#    journal_info3 = journalcheck(venue)[3]
    tmp = [[abs_words,number_authors,citation,number_ref,num_title,title,venue,venue1,year,authorimapct1,authorimapct2,authorpaper1,authorpaper2]] #,journal_info0,journal_info1,journal_info2,journal_info3]]
    df = pd.DataFrame(tmp,columns=['abs_words','number_authors','citation','number_ref','num_title',"title",'venue','venue1','year','authorimapct1','authorimapct2','authorpaper1','authorpaper2'])#'journal_info0','journal_info1','journal_info2','journal_info3'])
    if not os.path.isfile("D:/讨论班/report5/data/citationnumber0newupdate.csv"):
        df.to_csv("D:/讨论班/report5/data/citationnumber0newupdate.csv",header ='column_names')
    else: # else it exists so append without writing the header
        df.to_csv("D:/讨论班/report5/data/citationnumber0newupdate.csv",mode = 'a',header=False)
    print('do action'+str(i))
    
f = open("D:/讨论班/report5/data/citationnumber0newupdate.csv")
article_info = pd.read_csv(f)
f.close()

joint_info = pd.merge(article_info,sjrmerge,how='left',left_on='venue1',right_on='journaltitle1')
joint_info.to_csv("D:/讨论班/report5/data/citationjoint0.csv")    
#datafinal = datafinal.append(df)

#datafinal.to_csv("D:/讨论班/report5/data/citationnumber0.csv")

# Reads and converts json to dict.
