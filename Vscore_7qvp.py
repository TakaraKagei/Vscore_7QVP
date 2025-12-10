import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


#7qvpを読み込む
#177-253の行しか読み込んでいない
with open('/Users/mishimac9/Desktop/衝突可視化研究/PDB構造情報/7qvp.cif', 'r') as file:
    lines = file.readlines()[177:253]  
    
#タンパク質名を抽出
#polymer natの後ろに続く文字列を検索してリストに格納
#RNAは行指定で取り除いている
chain_name = re.findall(r'polymer\s+nat\s+\'([^\']*)\'', ''.join(lines))

#chain IDを抽出
# 477から1155行目までの行で、？を含んでいる行だけを抽出
with open('/Users/mishimac9/Desktop/衝突可視化研究/PDB構造情報/7qvp.cif', 'r', encoding='utf-8') as file:
    lines = file.readlines()[477:1168]
    chain_id_lines = [line.strip() for line in lines if '?' in line]

# 記号を排除し、,で区切られたふたつの文字列を別々の要素としてリストに格納
processed_chain_ids = []
for line in chain_id_lines:
    # 記号を排除
    cleaned_line = re.sub(r'[^\w,]', '', line)
    # ,がない場合は文字列の一文字目がLあるいはSで始まっている場合は,xxxを、そうでない場合はxxx,を追加
    if ',' not in cleaned_line:
        if cleaned_line.startswith(('L', 'S')):
            cleaned_line = cleaned_line + ',xxx'
        else:
            cleaned_line = 'xxx,' + cleaned_line
    # ,で区切られた文字列を分割してリストに追加
    split_line = cleaned_line.split(',')
    processed_chain_ids.extend(split_line)






#1111111111111111111111--------------------------------------------------------

#Entity IDとタンパク質名を対応づける(Key:Value = Entity ID:タンパク質名)
#Entity IDは9~84まで
#空の辞書に要素を追加していく
Entity_name = {}
x = 0
for x in range(len(chain_name)):
    Entity_name[x+9] = chain_name[x]

#Entity IDとchain IDを対応付する(Key:Value = Entity ID:chain ID)
#Entity IDは9~84まで
#空の辞書に要素を追加していく
#trailingとleadingは区別してchain IDを格納

#まずはtrailingとleadingをA, Bとして区別して、chainIDをそれぞれリストに格納
chainID_A = []
chainID_B = []

x = 0
for x in range(len(chain_id_lines)):
   chainID_A.append(processed_chain_ids[2*x])
   chainID_B.append(processed_chain_ids[2*x+1])

#Entity IDとchain IDの対応付けを行う(key:value = Entity ID:chain ID)
Entity_chainID_A = {}
Entity_chainID_B = {}

x = 0
for x in range(len(chain_id_lines)):
    Entity_chainID_A[x+9] = chainID_A[x]

x = 0
for x in range(len(chain_id_lines)):
    Entity_chainID_B[x+9] = chainID_B[x]



#2222222222222222222222--------------------------------------------------------
#座標を入れる空の辞書を作成
#N, C末端以外も一旦全部格納

#Entity ID 9~84までの,全原子座標が書かれた行を読み込む
with open('/Users/mishimac9/Desktop/衝突可視化研究/PDB構造情報/7qvp.cif', 'r') as file:
    lines = file.readlines()[289297:534222]

#chain IDで検索をかけ、ヒットした行を全てリストに格納
#アミノ酸の3文字表記とchain IDの文字列重複を避けるため、space + chain_ID + spaceを検索する
posi_A = {}
posi_B = {}

#通しアルファベットとchain IDの重複も避けるために、行に出てくるふたつ目の?以降の文字だけを検索対象としている
for n in range(len(chain_id_lines)):
    keyword_A = " " + chainID_A[n] + " "
    keyword_B = " " + chainID_B[n] + " "
    posi_A[chainID_A[n]] = [line for line in lines if keyword_A in line.split('?', 2)[-1]]
    posi_B[chainID_B[n]] = [line for line in lines if keyword_B in line.split('?', 2)[-1]]

#posi_A, posi_Bに格納した情報から、N末・C末のCAだけを抽出
#最初にヒットしたCAの行だけを抽出するようにする
#通しアルファベットにCAがあるため、"."以前の文字列のみを検索するようにする
keyword = "CA"
N_A = {}
C_A = {}
N_B = {}
C_B = {}

n = 0

keyword = "CA"
for n in range(len(chain_id_lines)):
    found = next((item for item in posi_A[chainID_A[n]] if keyword in item.split('.', 1)[0]), None)
    N_A[chainID_A[n]] = found
    found = next((item for item in reversed(posi_A[chainID_A[n]]) if keyword in item.split('.', 1)[0]), None)
    C_A[chainID_A[n]] = found
    found = next((item for item in posi_B[chainID_B[n]] if keyword in item.split('.', 1)[0]), None)
    N_B[chainID_B[n]] = found
    found = next((item for item in reversed(posi_B[chainID_B[n]]) if keyword in item.split('.', 1)[0]), None)
    C_B[chainID_B[n]] = found

#N_A, C_A, N_B, C_Bに格納された情報から、座標を含む部分だけを抽出
#座標の情報は、空白で区切られている
#?で囲まれた部分の情報だけを取り出す
#空白は","に置換
for n in range(len(chain_id_lines)):
    if N_A[chainID_A[n]]:
        N_A[chainID_A[n]] = re.findall(r'\?([^\?]*)\?', N_A[chainID_A[n]])[0].strip().replace(' ', ',')
    if C_A[chainID_A[n]]:
        C_A[chainID_A[n]] = re.findall(r'\?([^\?]*)\?', C_A[chainID_A[n]])[0].strip().replace(' ', ',')
    if N_B[chainID_B[n]]:
        N_B[chainID_B[n]] = re.findall(r'\?([^\?]*)\?', N_B[chainID_B[n]])[0].strip().replace(' ', ',')
    if C_B[chainID_B[n]]:
        C_B[chainID_B[n]] = re.findall(r'\?([^\?]*)\?', C_B[chainID_B[n]])[0].strip().replace(' ', ',')

#N_A, C_A, N_B, C_Bに格納された情報から、座標を抽出
#三つ目の,までを抽出
#3つ目の,は含まない
for n in range(len(chain_id_lines)):
    if N_A[chainID_A[n]]:
        N_A[chainID_A[n]] = N_A[chainID_A[n]].split(',', 3)[:3]
    if C_A[chainID_A[n]]:
        C_A[chainID_A[n]] = C_A[chainID_A[n]].split(',', 3)[:3]
    if N_B[chainID_B[n]]:
        N_B[chainID_B[n]] = N_B[chainID_B[n]].split(',', 3)[:3]
    if C_B[chainID_B[n]]:
        C_B[chainID_B[n]] = C_B[chainID_B[n]].split(',', 3)[:3]


#33333333333333333333333--------------------------------------------------------
#計算に入る
#先にDm(Distance in Monosome)を算出する
#AとBの両方でDmを求め、ふたつの値を平均した行列をNN, NC, CN, CCの4種類作成する

#AにおけるDmの算出
#NN, NC, CN, CCの4種類


m = 0
n = 0

Dm_NN_A = np.zeros((len(chainID_A), len(chainID_A)))

for m in range(len(chainID_A)):
    for n in range(len(chainID_A)):
        Dm_NN_A[m, n] = np.linalg.norm(np.array(N_A[chainID_A[m]], dtype=float) - np.array(N_A[chainID_A[n]], dtype=float))

Dm_NC_A = np.zeros((len(chainID_A), len(chainID_A)))

for m in range(len(chainID_A)):
    for n in range(len(chainID_A)):
        Dm_NC_A[m, n] = np.linalg.norm(np.array(N_A[chainID_A[m]], dtype=float) - np.array(C_A[chainID_A[n]], dtype=float))

Dm_CN_A = np.zeros((len(chainID_A), len(chainID_A)))

for m in range(len(chainID_A)):
    for n in range(len(chainID_A)):
        Dm_CN_A[m, n] = np.linalg.norm(np.array(C_A[chainID_A[m]], dtype=float) - np.array(N_A[chainID_A[n]], dtype=float))

Dm_CC_A = np.zeros((len(chainID_A), len(chainID_A)))

for m in range(len(chainID_A)):
    for n in range(len(chainID_A)):
        Dm_CC_A[m, n] = np.linalg.norm(np.array(C_A[chainID_A[m]], dtype=float) - np.array(C_A[chainID_A[n]], dtype=float))

#BにおけるDmの算出
#NN, NC, CN, CCの4種類

m = 0
n = 0

Dm_NN_B = np.zeros((len(chainID_B), len(chainID_B)))
                   
for m in range(len(chainID_B)):
    for n in range(len(chainID_B)):
        Dm_NN_B[m, n] = np.linalg.norm(np.array(N_B[chainID_B[m]], dtype=float) - np.array(N_B[chainID_B[n]], dtype=float))

Dm_NC_B = np.zeros((len(chainID_B), len(chainID_B)))
                   
for m in range(len(chainID_B)):
    for n in range(len(chainID_B)):
        Dm_NC_B[m, n] = np.linalg.norm(np.array(N_B[chainID_B[m]], dtype=float) - np.array(C_B[chainID_B[n]], dtype=float))

Dm_CN_B = np.zeros((len(chainID_B), len(chainID_B)))

for m in range(len(chainID_B)):
    for n in range(len(chainID_B)):
        Dm_CN_B[m, n] = np.linalg.norm(np.array(C_B[chainID_B[m]], dtype=float) - np.array(N_B[chainID_B[n]], dtype=float))

Dm_CC_B = np.zeros((len(chainID_B), len(chainID_B)))
                   
for m in range(len(chainID_B)):
    for n in range(len(chainID_B)):
        Dm_CC_B[m, n] = np.linalg.norm(np.array(C_B[chainID_B[m]], dtype=float) - np.array(C_B[chainID_B[n]], dtype=float))


#NN, NC, CN, CCの平均をとる
Dm_NN = np.nanmean([Dm_NN_A, Dm_NN_B], axis=0) 
Dm_NC = np.nanmean([Dm_NC_A, Dm_NC_B], axis=0) 
Dm_CN = np.nanmean([Dm_CN_A, Dm_CN_B], axis=0)
Dm_CC = np.nanmean([Dm_CC_A, Dm_CC_B], axis=0) 

# DataFrameを作成
df_Dm_NN = pd.DataFrame(Dm_NN, index=['N'] * len(chainID_A), columns=['N'] * len(chainID_A))
df_Dm_NC = pd.DataFrame(Dm_NC, index=['N'] * len(chainID_A), columns=['C'] * len(chainID_A))
df_Dm_CN = pd.DataFrame(Dm_CN, index=['C'] * len(chainID_A), columns=['N'] * len(chainID_A))
df_Dm_CC = pd.DataFrame(Dm_CC, index=['C'] * len(chainID_A), columns=['C'] * len(chainID_A))


#444444444444444444444444--------------------------------------------------------
#次にDc(Distance in Complex)を算出する
#NA_NB, NA_CB, CA_NB, CA_CBの4種類

Dc_NA_NB = np.zeros((len(chainID_A), len(chainID_B)))

for m in range(len(chainID_A)):
    for n in range(len(chainID_B)):
        Dc_NA_NB[m, n] = np.linalg.norm(np.array(N_A[chainID_A[m]], dtype=float) - np.array(N_B[chainID_B[n]], dtype=float))

Dc_NA_CB = np.zeros((len(chainID_A), len(chainID_B)))

for m in range(len(chainID_A)):
    for n in range(len(chainID_B)):
        Dc_NA_CB[m, n] = np.linalg.norm(np.array(N_A[chainID_A[m]], dtype=float) - np.array(C_B[chainID_B[n]], dtype=float))

Dc_CA_NB = np.zeros((len(chainID_A), len(chainID_B)))

for m in range(len(chainID_A)):
    for n in range(len(chainID_B)):
        Dc_CA_NB[m, n] = np.linalg.norm(np.array(C_A[chainID_A[m]], dtype=float) - np.array(N_B[chainID_B[n]], dtype=float))

Dc_CA_CB = np.zeros((len(chainID_A), len(chainID_B)))

for m in range(len(chainID_A)):
    for n in range(len(chainID_B)):
        Dc_CA_CB[m, n] = np.linalg.norm(np.array(C_A[chainID_A[m]], dtype=float) - np.array(C_B[chainID_B[n]], dtype=float))

# DataFrameを作成
df_Dc_NA_NB = pd.DataFrame(Dc_NA_NB, index=chainID_A, columns=chainID_B)
df_Dc_NA_CB = pd.DataFrame(Dc_NA_CB, index=chainID_A, columns=chainID_B)
df_Dc_CA_NB = pd.DataFrame(Dc_CA_NB, index=chainID_A, columns=chainID_B)
df_Dc_CA_CB = pd.DataFrame(Dc_CA_CB, index=chainID_A, columns=chainID_B)


#555555555555555555555555--------------------------------------------------------
#Dc, Dmを用いてVscoreを算出する
#Vscore = Dm / Dc

Vscore_NA_NB = Dm_NN / Dc_NA_NB
Vscore_NA_CB = Dm_NC / Dc_NA_CB
Vscore_CA_NB = Dm_CN / Dc_CA_NB
Vscore_CA_CB = Dm_CC / Dc_CA_CB


# DataFrameを作成
df_Vscore_NA_NB = pd.DataFrame(Vscore_NA_NB, index=chainID_A, columns=chainID_B)
df_Vscore_NA_CB = pd.DataFrame(Vscore_NA_CB, index=chainID_A, columns=chainID_B)
df_Vscore_CA_NB = pd.DataFrame(Vscore_CA_NB, index=chainID_A, columns=chainID_B)
df_Vscore_CA_CB = pd.DataFrame(Vscore_CA_CB, index=chainID_A, columns=chainID_B)

print(df_Vscore_NA_NB)
print(df_Vscore_NA_CB)
print(df_Vscore_CA_NB)
print(df_Vscore_CA_CB)


#666666666666666666666666--------------------------------------------------------
#順位表を作る
#Vscoreを降順に並べ替え、順位をつける
#計算に使ったchainIDを表示させる
# Vscoreを一列にまとめる
Vscore_combined = pd.concat([
    df_Vscore_NA_NB.stack().rename('Vscore').reset_index().assign(Type_A='N', Type_B='N'),
    df_Vscore_NA_CB.stack().rename('Vscore').reset_index().assign(Type_A='N', Type_B='C'),
    df_Vscore_CA_NB.stack().rename('Vscore').reset_index().assign(Type_A='C', Type_B='N'),
    df_Vscore_CA_CB.stack().rename('Vscore').reset_index().assign(Type_A='C', Type_B='C')
], axis=0)

# Vscoreを降順に並べ替え
Vscore_sorted = Vscore_combined.sort_values(by='Vscore', ascending=False)

# 列名を設定
Vscore_sorted.columns = ['chainID_A', 'chainID_B', 'Vscore', 'Type_A', 'Type_B']

# 順位をindexとして設定
Vscore_sorted.index = range(1, len(Vscore_sorted) + 1)

# 上位20位を表示
print(Vscore_sorted.head(40))

# Vscore が 3以上のものを抽出
filtered = Vscore_sorted[Vscore_sorted['Vscore'] >= 3].copy()

# 順位振り直し
filtered.reset_index(drop=True, inplace=True)
filtered.index += 1  # index: 1,2,3,...

# 色リスト
colors = []
for _, row in filtered.iterrows():
    a = row['chainID_A'][0]  # 先頭文字
    b = row['chainID_B'][0]

    if a in ['L','M'] and b in ['L','M']:
        colors.append("#98d98e")      # LSU / LSU (灰)
    elif a in ['S','R'] and b in ['S','R']:
        colors.append("#8AB6E6")      # SSU / SSU (青)
    else:
        colors.append("#E68A8A")      # mixed (緑)

plt.figure(figsize=(10,10))

# 横向き棒
plt.barh(filtered.index, filtered['Vscore'],
         color=colors, edgecolor='black', linewidth=1.0)

n = len(filtered)

# 縦軸（順位）の目盛り → 1,5,10,15,...
ticks = [1] + list(range(5, n+1, 5)) if n >= 5 else [1]
plt.yticks(ticks)

# 軸目盛り（数字）のフォントサイズ
plt.tick_params(axis='x', labelsize=12)  # 横軸の数字
plt.tick_params(axis='y', labelsize=12)  # 縦軸の数字

# 最小順位（1）が上に来るように反転
plt.gca().invert_yaxis()

# 枠線調整
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(1.5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_ylim(0.5, len(filtered) + 0.5)  # バーの上下の空白をなくす
ax.invert_yaxis()  # 順位1が上に来るように反転

# 凡例
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#E68A8A", label='LSU / SSU'),
    Patch(facecolor="#8AB6E6", label='SSU / SSU'),
    Patch(facecolor="#98d98e", label='LSU / LSU')
]
plt.legend(handles=legend_elements, loc='lower right', fontsize=15)

plt.xlabel("Vscore", fontsize=15)
plt.ylabel("Rank", fontsize=15)
plt.title("Vscore ≥ 3", fontsize=20)
plt.tight_layout()
plt.show()
