from tkinter import *
from tkinter import filedialog, messagebox
from tkinter.scrolledtext import ScrolledText
import time
import re

file_name = ''
codon_table = {
    'Ala':['GCG','GCA','GCT','GCC'],
    'Cys':['TGT','TGC'],
    'Asp':['GAT','GAC'],
    'Glu':['GAG' ,'GAA'],
    'Phe':['TTT' ,'TTC'],
    'Gly':['GGG','GGA','GGT','GGC'],
    'His':['CAT','CAC'],
    'Ile':['ATA','ATT','ATC'],
    'Lys':['AAG','AAA'],
    'Leu':['TTG','TTA','CTG','CTA','CTT','CTC'],
    'Met':['ATG'],
    'Asn':['AAT','AAC'],
    'Pro':['CCG','CCA','CCT','CCC'],
    'Gln':['CAG','CAA'],
    'Arg':['AGG','AGA','CGG','CGA','CGA','CGT','CGC'],
    'Ser':['AGT','AGC','TCG','TCA','TCT','TCC'],
    'Thr':['ACG','ACA','ACT','ACC'],
    'Val':['GTG','GTA','GTT','GTC'],
    'Trp':['TGG'],
    'Tyr':['TAT','TAC'],
    '*':['TGA','TAG','TAA']
}

# 读入fasta序列（代替SeqIO.parse）
def fasta_generator(file_path):
    with open(file_path, 'r') as file:
        first_line = file.readline().strip()
        if not first_line.startswith('>'):  # 给一个fasta格式判断
            raise ValueError('Invalid file format')
        sequence_id = first_line[1:]  # 储存当前序列名称
        sequence = ''   # 储存当前序列内容
        for line in file:
            line = line.strip()
            if line.startswith('>'):    # >开头为序列名称行
                if sequence_id:     # 如果已存在序列名称，返回上一个序列名称和内容
                    yield {'id': sequence_id, 'seq': sequence}
                sequence_id = line[1:]
                sequence = ''
            else:       # 否则，行为序列内容
                sequence += line
        if sequence_id:     # 处理完最后一个序列，返回最后一个序列名称和内容
            yield {'id': sequence_id, 'seq': sequence}

class MY_GUI():

    # 初始化函数，定义实例的属性
    def __init__(self,init_window_name):
        self.init_window_name = init_window_name
        self.image_file = None  # PhotoImage没有引用会自动销毁，这里需要显示引用

    # 定义一个类方法，设置窗口
    def set_init_window(self):
        self.init_window_name.title("密码子统计工具_v1.0")                    
        self.init_window_name.geometry('1048x680+400+150')    # 1068x681窗口大小，+横坐标 +纵坐标 定义窗口弹出时的默认展示位置
        self.init_window_name.attributes("-alpha",1)  # 虚化，值越小虚化程度越高
        self.init_window_name.iconbitmap('pic/bitbug_favicon.ico')
        self.init_window_name.resizable(0,0)    # 禁止改变窗口大小
        canvas = Canvas(self.init_window_name, width=1024, height=680, bg=None)
        self.image_file = PhotoImage(file="pic/bg.gif")
        self.resized_image = self.image_file.zoom(2, 2)
        canvas.create_image(520, 35, anchor='n', image=self.resized_image)
        canvas.grid(row=0,rowspan = 20,column=0,columnspan=23)
        # 标签
        self.init_data_label = Label(self.init_window_name, text="导入序列")
        self.init_data_label.grid(row=0, column=0)
        self.result_data_label = Label(self.init_window_name, text="Number(Frequency)")
        self.result_data_label.grid(row=0, column=12)
        self.result_out_label = Label(self.init_window_name, text="Codon Usage results")
        self.result_out_label.grid(row=6,column=13)
        self.log_label = Label(self.init_window_name, text="运行日志")
        self.log_label.grid(row=11, column=0)
        # 文本框
        self.init_data_Text = ScrolledText(self.init_window_name, width=60, height=35)  # 原始数据录入框
        self.init_data_Text.grid(row=1, column=0, rowspan=10, columnspan=10, padx=20,pady=5)
        self.init_data_Text.bind('<KeyPress>', lambda f: 'break')   # 绑定事件禁止键入
        self.log_data_Text = ScrolledText(self.init_window_name, width=60, height=9,)  # 日志框
        self.log_data_Text.grid(row=12, column=0, rowspan=5, columnspan=10,padx=20,pady=5)
        self.log_data_Text.bind('<KeyPress>', lambda f: 'break')
        self.result_data_Text = ScrolledText(self.init_window_name, width=60, height=20, wrap='none')  # 频数频率展示
        self.result_data_Text.grid(row=1, column=12, rowspan=5, columnspan=10, padx=10,pady=5)
        self.result_data_Text.bind('<KeyPress>', lambda f: 'break')
        self.result_fre_Text = ScrolledText(self.init_window_name, width=40, height=20)     # 结果展示
        self.result_fre_Text.grid(row=7, column=12, rowspan=10, columnspan=10,padx=10,pady=5)
        self.result_fre_Text.bind('<KeyPress>', lambda f: 'break')
        # 按钮
        self.load_file_button = Button(self.init_window_name, text = "导入fasta文件",borderwidth=2,relief=RAISED, command=self.select_file) # 绑定内部命令
        self.load_file_button.grid(row=2, column=11)
        self.codon_count_button = Button(self.init_window_name, text="开始统计",borderwidth=2,relief=RAISED,command=self.codon_count)
        self.codon_count_button.grid(row=3, column=11)
        self.export_count_button = Button(self.init_window_name, text="导出表格",borderwidth=2,relief=RAISED,command=self.export_count)
        self.export_count_button.grid(row=4, column=11)
        self.export_fre_button = Button(self.init_window_name, text="导出结果", borderwidth=2,relief=RAISED,command=self.export_fre)
        self.export_fre_button.grid(row=5,column=11)
        self.clean_button = Button(self.init_window_name, text='清空窗口', borderwidth=2,relief=RAISED,command=self.window_clean)
        self.clean_button.grid(row=6,column=11)
        # 滚轮与绑定（ScrolledText只带有垂直滚动条）
        self.scrollbar_x = Scrollbar(self.init_window_name, orient=HORIZONTAL)  # 创建滚动条部件
        self.result_data_Text.config(xscrollcommand=self.scrollbar_x.set)   # 文本框-控制-滚动条
        self.scrollbar_x.config(command=self.result_data_Text.xview)    # 滚动条-控制-文本框
        self.scrollbar_x.grid(row=5, column=12, rowspan=1,columnspan=10, sticky="ESW",padx=10)   # 设置滚动条位置

    # 导入序列功能函数
    def select_file(self):
        global file_name
        file_name = filedialog.askopenfilename(title="选择文件")
        if file_name != '':
            try:  
                records = fasta_generator(file_name)
                number = 0  # 计数，导入多少序列
                for i in records:
                    name = i['id']
                    seq = i['seq']
                    if len(seq)>48:
                        seq = seq[:45] + '...' + seq[-3:]
                    self.init_data_Text.insert(END,f"Name:{name}\nSeq:'{seq}'\n")
                    number += 1
                self.write_log_to_Text(f"INFO:fasta文件导入成功！共加载{number}条序列。")
            except:
                messagebox.showerror("ERROR", "载入失败，请检查文件是否为fasta格式！")
                self.write_log_to_Text("ERROR:载入失败，请检查文件格式。")

    # 分析统计功能函数
    def codon_count(self):
        if file_name != '':
            records = fasta_generator(file_name)
            CodonsDict = {codon: 0 for codon_list in codon_table.values() for codon in codon_list}  # 新的字典，统计密码子数量用
            # 输出表头
            self.result_data_Text.insert(END, f"{'Name':<15}")
            for key in CodonsDict:
                    self.result_data_Text.insert(END, f'{key:<12}')
            self.result_data_Text.insert(END,'\n')
            for i in records:
                # 每条序列判断ATG开头，是否有屏蔽序列,长度是否为3的倍数
                if i['seq'].startswith('ATG') and 'N' not in i['seq'] and len(i['seq']) % 3 ==0:
                    for j in range(0, len(str(i['seq'])), 3):
                        codon = str(i['seq'])[j:j+3]
                        if codon in CodonsDict.keys():
                            CodonsDict[codon] +=1
                        else:
                            self.write_log_to_Text("WARNING:序列%s存在未识别的密码子，跳过。" % (i['id']))
                            break
                    total = sum([CodonsDict[key] for key in CodonsDict.keys()])
                    name = i['id']
                    self.result_data_Text.insert(END, f'{name:<15}')    # 这里的f-string格式化输出f'{i['id']:<15}'会报错，所以用了个变量name代替
                    self.result_fre_Text.insert(END,'Results for %d residue sequence "%s":\n\nAA\tCodon\tNumber\tFrequency\n\n' % (total, name))
                    # 计算频率
                    for key, value in CodonsDict.items():
                        frequency = '%.2f' % (value * 3000 / total)
                        content = '%d(%s)' % (value, frequency)
                        self.result_data_Text.insert(END,f'{content:<12}')
                        # 结果文本框内的输出
                        for key_, value_ in codon_table.items():
                            if key in value_:
                                AA = key_
                                self.result_fre_Text.insert(END, '%s\t%s\t%d\t%s\n' % (AA, key, value, frequency))
                    self.result_fre_Text.insert(END,'\n----------------------------------------\n' )
                    self.result_data_Text.insert(END,'\n')
                else:
                    self.write_log_to_Text("WARNING:序列%s非ATG开头/存在屏蔽序列/非3的倍数，该序列将不会出现在统计结果中。" % (i['id']))
            self.write_log_to_Text("INFO:统计结束!")
        else:
            self.write_log_to_Text("ERROR:请先载入fasta文件！")

    # 导出频数频率统计表
    def export_count(self):
        content = self.result_data_Text.get('1.0', END)
        if 'Name' not in content:
            self.write_log_to_Text("ERROR:请先点击“载入fasta文件”并点击“开始统计”按钮！")
        else:
            try:
                count_file = filedialog.asksaveasfilename()
                count_file = count_file + '.csv'
                with open(count_file, 'w') as count:
                    content = re.sub(r"[^\S\r\n]+",',',content) # 正则匹配换行符之外的所有空格，处理成csv格式的输出
                    count.write(content)
                self.write_log_to_Text("INFO:密码子频数频率统计表保存成功!文件路径：%s" % (count_file))
            except:
                self.write_log_to_Text("ERROR:ERROR:导出失败，请检查是否有同名文件未关闭！")

    # 导出结果文件
    def export_fre(self):
        content = self.result_fre_Text.get('1.0', END)
        if 'Results' not in content:
            self.write_log_to_Text("ERROR:请先点击“载入fasta文件”并点击“开始统计”按钮！")
        else:
            try:
                result_file = filedialog.asksaveasfilename()
                result_file = result_file + '.txt'
                with open(result_file, 'w') as out:
                    out.write(content)
                self.write_log_to_Text("INFO:结果文件保存成功!保存路径：%s" % (result_file))
            except:
                self.write_log_to_Text("ERROR:导出失败，请检查是否有同名文件未关闭！")

    # 清空所有文本框
    def window_clean(self):
        self.init_data_Text.delete(1.0, END)
        self.log_data_Text.delete(1.0, END)
        self.result_data_Text.delete(1.0, END)
        self.result_fre_Text.delete(1.0, END)

    # 日志打印
    def write_log_to_Text(self,logmsg):
        current_time = time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time()))
        logmsg_in = str(current_time) +" " + str(logmsg) + "\n"
        self.log_data_Text.insert(END, logmsg_in)

def gui_start():
    init_window = Tk()
    GUI = MY_GUI(init_window)
    GUI.set_init_window()
    init_window.mainloop()

gui_start()