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

# GUIのためにPySimpleGUIモジュールをインポートする
import PySimpleGUI as sg

# ラジオボタン周りのフレームの設定
frame1 = [  [sg.Radio('SpCas9 + pSNR52-sgRNA (plasmid 15-13)', 1, key='-SpCas9_pSNR52-', default=True)],
            [sg.Radio('SpCas9 + pGAL1-sgRNA (plasmid 16-15)', 1, key='-SpCas9_pGAL1-')],
            [sg.Radio('SaCas9 + pGAL1-sgRNA (plasmid 17-31)', 1, key='-SaCas9-')],
            [sg.Radio('enAsCas12a + pGAL1-crRNA (plasmid 16-16)', 1, key='-enAsCas12a-')]]

# レイアウト
layout = [  [sg.Text('Input target name(s) and target sequence(s):')],
            [sg.Multiline(size=(100, 10), key='textbox1')],
            [sg.Button('Clear the input')],
            [sg.Frame('Cas protein + sgRNA/crRNA', frame1)],
            [sg.Button('Design oligo DNA sequences')],
            [sg.Text('Oligo DNA sequences for Golden Gate Assembly are shown below:')],
            [sg.Text('(You can copy the result to a spread sheet application such as Microsoft Excel.)')],
            [sg.Multiline(size=(100, 10), key='textbox2')],
            [sg.Button('Clear the result')]]

# ウィンドウを作成する
window = sg.Window('Oligo DNA designer for budding yeast genome-editing plasmid construction (v.200916)', layout).Finalize()

# イベントループ
while True:
    event, values = window.read()
    HHv = ''
    HHc = ''

    # 'Design oligo DNA sequences'ボタンが押されたときの動作
    if event in ('Design oligo DNA sequences'):

        # ラジオボタン'SpCas9 + pSNR52-sgRNA (plasmid 15-13)'が押されていたときの動作
        if values['-SpCas9_pSNR52-'] == True:

            # 変数 HHc
            # hammerhead ribozymeの後半の37塩基（固定配列）
            HHc = 'CTGATGAGTCCGTGAGGACGAAACGAGTAAGCTCGTC'

            # 変数 GGA5
            # Golden Gate Assembly向けにFwdオリゴの5'末に付加する配列
            GGA5 = 'GATC'

            # 変数 GGA3
            # Golden Gate Assembly向けにRevオリゴの5'末に付加する配列
            GGA3 = 'AAAC'

            # 変数 Header
            # 最初の行に出力する文字列
            Header = 'Target name\tTarget seq\tHH + Target seq\tFwd seq for GGA (15-13)\tRev seq for GGA (15-13)'
            # 変数 Target_name
            # ターゲット配列名
            Target_name = ''

            # 変数 Target
            # ターゲット配列
            Target = ''

            # 入力テキストボックス内のテキストを読み出す
            input_text = values['textbox1']

            # 出力テキストボックスにヘッダーを書き込む
            window['textbox2'].print(Header)

            # 入力テキストボックス内のテキストを1行ごとに分割し、オリゴDNA配列をデザインする
            lines = input_text.split('\n')
            for line in lines:
                line_strip = line.strip()
                line_split = line_strip.split()
                if len(line_split) == 2:
                    Target_name = line_split[0]
                    Target = line_split[1]
                    # 変数 HHv
                    # hammerhead ribozymeの最初の6塩基（後続の配列に依存して変化する）
                    HHv = HHvariable(Target)
                    HH = HHv + HHc
                    HH_Target = HH + Target
                    F_seq_GGA = GGA5 + HH_Target
                    R_seq_GGA = GGA3 + complement_seq(HH_Target)
                    window['textbox2'].print(Target_name + '\t'+ Target +'\t' + HH_Target + '\t' + F_seq_GGA + '\t' + R_seq_GGA)

        # ラジオボタン'SpCas9 + pGAL1-sgRNA (plasmid 16-15)'が押されていたときの動作
        if values['-SpCas9_pGAL1-'] == True:

            HHc = 'CTGATGAGTCCGTGAGGACGAAACGAGTAAGCTCGTC'
            GGA5 = 'GGAG'
            GGA3 = 'AAAC'
            Header = 'Target name\tTarget seq\tHH + Target seq\tFwd seq for GGA (16-15)\tRev seq for GGA (16-15)'
            # 変数 Target_name
            # ターゲット配列名
            Target_name = ''

            # 変数 Target
            # ターゲット配列
            Target = ''

            # 入力テキストボックス内のテキストを読み出す
            input_text = values['textbox1']

            # 出力テキストボックスにヘッダーを書き込む
            window['textbox2'].print(Header)

            # 入力テキストボックス内のテキストを1行ごとに分割し、オリゴDNA配列をデザインする
            lines = input_text.split('\n')
            for line in lines:
                line_strip = line.strip()
                line_split = line_strip.split()
                if len(line_split) == 2:
                    Target_name = line_split[0]
                    Target = line_split[1]
                    HHv = HHvariable(Target)
                    HH = HHv + HHc
                    HH_Target = HH + Target
                    F_seq_GGA = GGA5 + HH_Target
                    R_seq_GGA = GGA3 + complement_seq(HH_Target)
                    window['textbox2'].print(Target_name + '\t'+ Target +'\t' + HH_Target + '\t' + F_seq_GGA + '\t' + R_seq_GGA)

        # ラジオボタン'SaCas9 + pGAL1-sgRNA (plasmid 17-31)'が押されていたときの動作
        if values['-SaCas9-'] == True:
            HHc = 'CTGATGAGTCCGTGAGGACGAAACGAGTAAGCTCGTC'
            GGA5 = 'GGAG'
            GGA3 = 'TAAC'
            Header = 'Target name\tTarget seq\tHH + Target seq\tFwd seq for GGA (17-31)\tRev seq for GGA (17-31)'
            # 変数 Target_name
            # ターゲット配列名
            Target_name = ''

            # 変数 Target
            # ターゲット配列
            Target = ''

            # 入力テキストボックス内のテキストを読み出す
            input_text = values['textbox1']

            # 出力テキストボックスにヘッダーを書き込む
            window['textbox2'].print(Header)

            # 入力テキストボックス内のテキストを1行ごとに分割し、オリゴDNA配列をデザインする
            lines = input_text.split('\n')
            for line in lines:
                line_strip = line.strip()
                line_split = line_strip.split()
                if len(line_split) == 2:
                    Target_name = line_split[0]
                    Target = line_split[1]
                    HHv = HHvariable(Target)
                    HH = HHv + HHc
                    HH_Target = HH + Target
                    F_seq_GGA = GGA5 + HH_Target
                    R_seq_GGA = GGA3 + complement_seq(HH_Target)
                    window['textbox2'].print(Target_name + '\t'+ Target +'\t' + HH_Target + '\t' + F_seq_GGA + '\t' + R_seq_GGA)

        # ラジオボタン'enAsCas12a + pGAL1-crRNA (plasmid 16-16)'が押されていたときの動作
        if values['-enAsCas12a-'] ==True:
            HHv = ''
            HHc = ''
            GGA5 = 'AGAT'
            GGA3 = 'AAAA'
            Header = 'Target name\tTarget seq\tHH + Target seq\tFwd seq for GGA (16-16)\tRev seq for GGA (16-16)'
            # 変数 Target_name
            # ターゲット配列名
            Target_name = ''

            # 変数 Target
            # ターゲット配列
            Target = ''

            # 入力テキストボックス内のテキストを読み出す
            input_text = values['textbox1']

            # 出力テキストボックスにヘッダーを書き込む
            window['textbox2'].print(Header)

            # 入力テキストボックス内のテキストを1行ごとに分割し、オリゴDNA配列をデザインする
            lines = input_text.split('\n')
            for line in lines:
                line_strip = line.strip()
                line_split = line_strip.split()
                if len(line_split) == 2:
                    Target_name = line_split[0]
                    Target = line_split[1]
                    HH = HHv + HHc
                    HH_Target = HH + Target
                    F_seq_GGA = GGA5 + HH_Target
                    R_seq_GGA = GGA3 + complement_seq(HH_Target)
                    window['textbox2'].print(Target_name + '\t'+ Target +'\t' + HH_Target + '\t' + F_seq_GGA + '\t' + R_seq_GGA)

    # 'Clear the input'ボタンが押されたときの動作
    if event in ('Clear the input'):
        values['textbox1'] = ''
        window['textbox1'].update(values['textbox1'])

    # 'Clear the result'ボタンが押されたときの動作
    if event in ('Clear the result'):
        values['textbox2'] = ''
        window['textbox2'].update(values['textbox2'])

    # ウィンドウ右上の×ボタンが押されたときの動作
    if event in (None, 'Close Window'):
        break

window.close()
