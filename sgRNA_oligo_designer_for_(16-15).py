# 関数 complement_nt
# 1塩基 x を受け取ってその相補的な塩基を返す
def complement_nt(x):
    if x == 'A':
        return 'T'
    elif x == 'a':
        return 'T'
    elif x == 'C':
        return 'G'
    elif x == 'c':
        return 'G'
    elif x == 'G':
        return 'C'
    elif x == 'g':
        return 'C'
    elif x == 'T':
        return 'A'
    elif x == 't':
        return 'A'
    else:
        return 'N'

# 関数 complement_seq
# 塩基配列 seq を受け取ってその相補的な配列を返す
def complement_seq(seq):
    n = len(seq)
    complement_seq =''
    while n>0:
        complement_seq = complement_seq + complement_nt(seq[n-1])
        n = n-1
    return complement_seq

# 関数 HHvariable
# 塩基配列 target を受け取ってその最初の6塩基に対するhammerhead ribozymeの最初の6塩基を返す
def HHvariable(target):
    HHvariable =''
    for n in range(6):
        HHvariable = complement_nt(target[n]) + HHvariable
    return HHvariable

# 変数 HHv
# hammerhead ribozymeの最初の6塩基（後続の配列に依存して変化する）
HHv = ''

# 変数 HHc
# hammerhead ribozymeの後半の37塩基（固定配列）
HHc = 'CTGATGAGTCCGTGAGGACGAAACGAGTAAGCTCGTC'

# 変数 GGA5
# Golden Gate Assembly向けにFwdオリゴの5'末に付加する配列
GGA5 = 'GGAG'

# 変数 GGA3
# Golden Gate Assembly向けにRevオリゴの5'末に付加する配列
GGA3 = 'AAAC'

# 変数 Header
# 最初の行に出力する文字列
Header = 'Target name\tTarget seq\tHH + Target seq\tFwd seq for GGA (16-15)\tRev seq for GGA (16-15)'

# argvを取得するためにsysモジュールをインポートする
import sys
# コマンドライン引数をargvs（リスト）に格納する
argvs = sys.argv
# コマンドライン引数の数を変数argcに格納する
argc = len(argvs)

# 入力するファイル群が指定されていないときは使い方を表示して終了する
if argc < 2:
	print("Usage: python {} [input_file.tsv] ".format(argvs[0]))
else:

    # 変数 File_name
    # ターゲット配列名とターゲット配列を記載したファイル（タブ区切り）のファイル名
    File_name = argvs[1]

    # 変数 File_name_output
    # 出力ファイル名（入力ファイル名の末尾に'.output.txt'を付加したもの）
    File_name_output = File_name + '.output.txt'

    # 変数 Target_name
    # ターゲット配列名
    Target_name = ''

    # 変数 Target
    # ターゲット配列
    Target = ''

    # 出力ファイルの先頭にヘッダーを書き込む
    with open(File_name_output, mode='a') as output:
        output.write(Header)
        output.write('\n')
        output.close()

    # 出力ファイルの2行目以降に配列情報を書き込む
    with open(File_name) as f:
        lines = f.readlines()
        for line in lines:
            line_strip = line.strip()
            line_split = line_strip.split()
            Target_name = line_split[0]
            Target = line_split[1]
            HHv = HHvariable(Target)
            HH = HHv + HHc
            HH_Target = HH + Target
            F_seq_GGA = GGA5 + HH_Target
            R_seq_GGA = GGA3 + complement_seq(HH_Target)
            with open(File_name_output, mode='a') as output:
                output.write(Target_name + '\t'+ Target +'\t' + HH_Target + '\t' + F_seq_GGA + '\t' + R_seq_GGA)
                output.write('\n')
                output.close()
            f.close()

    print('Output saved as ' + File_name_output)
