'''______________________________________________________________
txtファイルからfloat型の行列を作成します。
必要があれば下のcsv_convで行列へ変換します。

csv_conv関数で任意のファイルからリストを作成します。
クォーテーションは消えて数値の行列が作成されます。

_________________________________________________________________'''

import csv
import matplotlib.pyplot as plt
import math
import numpy as np
import random
from mpl_toolkits.mplot3d import Axes3D

def txt_csv_conv(filename): #入力はstr型で行う　拡張子を含まないように（含んでもできる）
    putname = str(filename).replace(".txt", "") + ".txt"
    output_list = []
    tot = []
    with open(putname,"r") as conv_file:
        mk = csv.reader(conv_file)
        header = next(mk) #headerを無視
        ne = list(mk)
        #del(ne[0]) #不要なデータを消去(headerのこと)
    for j in range(len(ne)):
        for  k in range(int(len(ne[j]))):
            a = float(ne[j][k])
            tot.append(a)
        output_list.append(tot)
        tot = []         #ここで文字列型をfloat型に変換
    return output_list

def p_selecter(filename,n):#nは粒子種を示す番号
    select_pmtx=[]
    p=txt_csv_conv(filename)
    if n==1:
        select_pmtx=p[0]
    else:
        select_pmtx=p[1]
    return select_pmtx

def pn_selecter(filename,n):
    select_pmtx=[]
    p=txt_csv_conv(filename)
    if n==1:
        select_pmtx=p[-2]
    else:
        select_pmtx=p[-1]
    return select_pmtx

"""=====CSV converter===========================================================
def csv_conv(file_name):  #入力はStr型で行ってください 拡張子を含んでもできるが含まないようにしてください
    putname = str(file_name).replace(".csv", "") + ".csv"
    output_list = []
    tot = []
    with open(putname,"r") as conv_file:
        mk = csv.reader(conv_file)
        ne = list(mk)
        del(ne[0])
       # ne.reverse()     #ここまででCSVファイルを行列に変換
    for j in range(len(ne)):
        for  k in range(int(len(ne[j]))):
            a = float(ne[j][k])
            tot.append(a)
        output_list.append(tot)
        tot = []         #ここで文字列型をfloat型に変換
    return output_list
==============================================================================="""

D_rhme = []
D_rhmp = []
D_ne = []
D_np = []
D_oe = []
D_op = []
D_pe = []
D_pp = []
Drhmrhm=[]

###########値の代入################################
D_rhme=p_selecter("/Users/anzairyoukei/git/task/fp.ota/dat/Drhmrhm.txt",1)
D_rhmp=p_selecter("/Users/anzairyoukei/git/task/fp.ota/dat/Drhmrhm.txt",2)

D_ne=p_selecter("/Users/anzairyoukei/git/task/fp.ota/dat/Dnewba.txt",1)
D_np=p_selecter("/Users/anzairyoukei/git/task/fp.ota/dat/Dnewba.txt",2)

D_oe=p_selecter("/Users/anzairyoukei/git/task/fp.ota/dat/Drwav.txt",1)
D_op=p_selecter("/Users/anzairyoukei/git/task/fp.ota/dat/Drwav.txt",2)

D_pe=p_selecter("/Users/anzairyoukei/git/task/fp.ota/dat/Dnewpla.txt",1)
D_pp=p_selecter("/Users/anzairyoukei/git/task/fp.ota/dat/Dnewpla.txt",2)
r = p_selecter("/Users/anzairyoukei/git/task/fp.ota/dat/rm.txt",1)
###########描画###################################
plt.plot(r, D_oe,label="D_ota")
plt.plot(r, D_ne,label="Neoclassical_ba")
#plt.plot(r, D_pe,label="Neoclassical_pla")
plt.xlabel('r/a', fontname="Times New Roman")
plt.ylabel('$D_e$[$m^2$/s]',fontname="Times New Roman")
plt.grid()
plt.legend()
plt.savefig("/Users/anzairyoukei/Desktop/径方向拡散_電子.png",dpi=100, pad_inches = 'tight')
plt.show()

plt.plot(r, D_op,label="D_ota")
plt.plot(r, D_np,label="Neoclassical_ba")
#plt.plot(r, D_ne,label="Neoclassical_pla")
plt.xlabel('r/a',fontname="Times New Roman")
plt.ylabel('$D_p$[$m^2$/s]',fontname="Times New Roman")
plt.grid()
plt.legend()
plt.savefig("/Users/anzairyoukei/Desktop/径方向拡散_陽子.png",dpi=100, pad_inches = 'tight')
plt.show()

print("=============================================================")
plt.plot(r, D_rhme,label="FOW")
plt.plot(r, D_oe,label="D_ota")
plt.plot(r, D_ne,label="Neoclassical_ba")
#plt.plot(r, D_pe,label="Neoclassical_pla")
plt.xlabel('r/a', fontname="Times New Roman")
plt.ylabel('$D_e$[$m^2$/s]',fontname="Times New Roman")
plt.ylim(0, 0.001)
plt.grid()
plt.legend()
#画像の保存と保存先について
plt.savefig("/Users/anzairyoukei/Desktop/径方向拡散_電子1.png",dpi=100, pad_inches = 'tight')
plt.show()

plt.plot(r, D_rhmp,label="FOW")
plt.plot(r, D_op,label="D_ota")
plt.plot(r, D_np,label="Neoclassical_ba")
#plt.plot(r, D_pp,label="Neoclassical_pla")
plt.xlabel('r/a',fontname="Times New Roman")
plt.ylabel('$D_p$[$m^2$/s]',fontname="Times New Roman")
plt.ylim(0, 0.001)
plt.grid()
plt.legend()
#画像の保存と保存先について
plt.savefig("/Users/anzairyoukei/Desktop/径方向拡散_陽子1.png",dpi=100, pad_inches = 'tight')
plt.show()
