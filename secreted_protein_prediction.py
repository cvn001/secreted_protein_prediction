import os
import sys
import time
import shutil
# import requests
import subprocess
import traceback
import matplotlib
import logging.handlers
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib import rcParams
from argparse import ArgumentParser
from collections import defaultdict, OrderedDict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool


def parse_cmdline():
    """
    Parse command-line arguments for script.
    :return: Input command-line arguments
    """
    parser = ArgumentParser(prog="secreted_protein_prediction.py")
    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        action="store",
        required=True,
        help="Input an annotation combine xlsx format file [Annotation_combine.xlsx]."
    )
    parser.add_argument(
        "-s",
        "--seq",
        dest="sequence",
        required=True,
        help="Input a reference protein sequence file [Homo_sapiens_9606_SP_20200509.fasta]."
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        action="store",
        default='Results',
        help="Output directory [default: Results]."
    )
    parser.add_argument(
        "-t",
        "--taxon",
        dest="taxonomy",
        action="store",
        required=True,
        choices=['archaea', 'gram+', 'gram-', 'animal', 'plant'],
        help="Organism [Archaea: 'archaea', Gram-positive: 'gram+', "
             "Gram-negative: 'gram-', Animal: 'animal' or Plant: 'plant']."
    )
    parser.add_argument(
        "-c",
        "--cpu",
        type=int,
        action="store",
        dest="threads",
        default=50,
        help="How many CPU threads will be used? No more than 80 [default 50]."
    )
    parser.add_argument(
        "-nu",
        "--no_uniprot",
        dest="no_uniprot",
        action="store_true",
        default=False,
        help="Whether to use Uniprot keywords to do prediction [default True]."
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        action="store_true",
        default=False,
        help="Give verbose output."
    )
    return parser.parse_args()


def last_exception():
    """ 返回上一个错误信息，用于logging打印出来"""
    exc_type, exc_value, exc_traceback = sys.exc_info()
    return ''.join(
        traceback.format_exception(exc_type, exc_value, exc_traceback))


def new_logger():
    """
    本函数用于生成屏幕日志
    """
    # Set up logging
    new_log = logging.getLogger('Start secreted proteins prediction: %s' % time.asctime())
    # 指定logging输出的格式
    # new_formatter = logging.Formatter('%(levelname)s: %(message)s')
    new_format = logging.Formatter('[%(levelname)s] - %(message)s')
    # 指定日志的最低输出级别
    new_log.setLevel(logging.DEBUG)
    # 控制台日志
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.formatter = new_format
    if args.verbose:
        console_handler.setLevel(logging.INFO)
    else:
        console_handler.setLevel(logging.WARNING)
    new_log.addHandler(console_handler)
    return new_log, new_format


def secretomep_runner(program, files):
    seq_file = files[0]
    result_file = files[1]
    cmd = '{} -s {} > {}'.format(program, seq_file, result_file)
    devnull = open(os.devnull, 'w')
    subprocess.call(cmd, shell=True, stdout=devnull, stderr=devnull)


def multiple_secretome(program, input_list, threads):
    pool = Pool(int(threads))
    for files in input_list:
        pool.apply_async(secretomep_runner, args=(program, files))
    pool.close()
    pool.join()


class MainPipeline:
    """
    Secreted protein analysis
    """

    def __init__(self):
        self.annotation_file = args.input
        self.df = None
        self.threads = str(args.threads)
        self.no_uniprot = args.no_uniprot
        self.temp_dir = os.path.join(output_dir, 'temp')
        os.makedirs(self.temp_dir, exist_ok=True)
        self.organism = args.taxonomy
        self.protein_dict = defaultdict(list)
        self.annotation_dict = defaultdict(list)
        self.classical_dict = defaultdict(list)
        self.non_classical_dict = defaultdict(list)
        self.reference_seq_file = args.sequence
        self.protein_seq_file = os.path.join(self.temp_dir, 'proteins.fasta')
        self.signalp = 'signalp'
        self.secretomep = '/home/lxc/bin/secretomep-1.0/secretomep'
        self.temp_result_file = ''
        self.excel_file = os.path.join(output_dir, 'Secreted_proteins_prediction.xlsx')
        self.statistics_dict = defaultdict()
        self.sub_file_list = []

    def load_data(self):
        """
        此函数用于导入Annotation_combine.xlsx中的蛋白和GO注释数据
        """
        logger.info('Loading input data.')
        tmp_df = pd.read_excel(self.annotation_file)
        tmp_df = tmp_df.replace(np.nan, '', regex=True)
        if self.no_uniprot:
            self.df = tmp_df[['Protein accession', 'Cellular Component']]
        else:
            if 'Uniprot Keywords' in tmp_df.columns:
                self.df = tmp_df[['Protein accession', 'Cellular Component', 'Uniprot Keywords']]
            else:
                logger.error("No [Uniprot Keywords] column found in Annotation_combine.xlsx! Aborting.")
                sys.exit(1)
        for index, row in tmp_df.iterrows():
            row_protein = row['Protein accession']
            row_go = row['Cellular Component'].strip()
            if self.no_uniprot:
                keywords = ''
            else:
                keywords = row['Uniprot Keywords'].strip()
            row_description = row['Protein description'].replace('"', '')
            row_name = row['Gene name']
            self.annotation_dict[row_protein] = [row_protein, row_description, row_name]
            if row_protein not in self.protein_dict:
                self.protein_dict[row_protein].append(row_go)
                self.protein_dict[row_protein].append(keywords)
        protein_number = len(self.protein_dict)
        logger.info('Total protein number: [{}].'.format(protein_number))

    def fetch_protein_sequences(self):
        """
        此函数用于提取全部鉴定蛋白的序列
        """
        logger.info('Fetching identified protein sequences.')
        records = []
        for seq_record in SeqIO.parse(self.reference_seq_file, 'fasta'):
            seq_id = str(seq_record.id)
            if seq_id in self.protein_dict:
                new_record = SeqRecord(id=seq_id, seq=seq_record.seq, description='')
                records.append(new_record)
        SeqIO.write(sequences=records, handle=self.protein_seq_file, format='fasta')

    def seqs_separation(self):
        """
        此函数按照进程数将原始的protein.fasta文件进行等比拆分，
        得到相同数目的子序列文件，便于并行运算。
        :return: parameters为一个列表(list)，包含三个元素，
        """
        part = int(self.threads)
        seq_dict = OrderedDict()
        i = 0
        for record in SeqIO.parse(self.protein_seq_file, 'fasta'):
            seq_dict[i] = SeqRecord(seq=record.seq, id=record.id, description='')
            i += 1
        div_groups = np.array_split(list(seq_dict.keys()), part)
        j = 0
        for group in div_groups:
            if len(group) > 0:
                group_records = []
                group_seq_file = os.path.join(self.temp_dir, 'sub_{}.fasta'.format(j))
                group_result_file = os.path.join(self.temp_dir, 'secretomep_sub_{}.txt'.format(j))
                self.sub_file_list.append([group_seq_file, group_result_file])
                for seq_index in group:
                    group_records.append(seq_dict[seq_index])
                SeqIO.write(sequences=group_records, handle=group_seq_file, format='fasta')
                j += 1

    # def uniprot_request(self):
    #     """
    #     此函数用于批量爬取Uniprot Keywords数据，该过程比较耗时
    #     """
    #     if self.no_uniprot:
    #         keywords = ''
    #         for each_protein in self.protein_dict.keys():
    #             self.protein_dict[each_protein].append(keywords)
    #     else:
    #         logger.info('Requesting Uniprot keywords, may take a while.')
    #         base = 'http://www.uniprot.org'
    #         kb_endpoint = '/uniprot/'
    #         web_headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; WOW64) '
    #                        'AppleWebKit/537.36 (KHTML, like Gecko) '
    #                        'Chrome/55.0.2883.87 Safari/537.36'}
    #         for each_protein in self.protein_dict.keys():
    #             payload = {'query': each_protein, 'format': 'tab', 'columns': 'keywords'}
    #             result = requests.get(base + kb_endpoint, params=payload, headers=web_headers)
    #             if result.ok:
    #                 result_line = result.text
    #                 k_list = result_line.strip().split('\n')
    #                 if len(k_list) > 1:
    #                     keywords = k_list[1]
    #                 else:
    #                     keywords = ''
    #             else:
    #                 keywords = ''
    #                 logger.warning('[{}] Uniprot keywords downloading failed.'.format(each_protein))
    #             self.protein_dict[each_protein].append(keywords)
    #             # 设置0.1秒钟的延迟对于Uniprot反爬虫机制更加安全
    #             time.sleep(0.1)

    def run_signalp(self):
        """
        此函数用于运行SignalP v5.0软件进行信号肽预测
        """
        signalp_file_prefix = os.path.join(self.temp_dir, 'result')
        if self.organism in ['plant', 'animal']:
            organism = 'euk'
        else:
            organism = self.organism
        cmd = '{} -fasta {} -org {} -format short -plot none -prefix {} ' \
              '-verbose=false'.format(self.signalp, self.protein_seq_file,
                                      organism, signalp_file_prefix)
        try:
            logger.info("Running SignalP v5.0b to predict classical signal peptides.")
            subprocess.call(cmd, shell=True)
        except OSError:
            logger.error("Try to run SignalP but failed! Aborting.")
            sys.exit(1)
        result_file = os.path.join(self.temp_dir, 'result_summary.signalp5')
        r_dict = defaultdict()
        with open(result_file, 'r') as f:
            for each_line in f.readlines():
                if '#' not in each_line:
                    r_list = each_line.strip('\r|\n').split('\t')
                    protein = r_list[0]
                    prediction = r_list[1]
                    r_dict[protein] = prediction
        for each_protein in self.protein_dict.keys():
            if each_protein in r_dict:
                each_prediction = r_dict[each_protein]
            else:
                each_prediction = 'OTHER'
            self.protein_dict[each_protein].append(each_prediction)

    def run_secretomep(self):
        """
        此函数用于运行SecretomeP软件预测非经典分泌蛋白
        """
        if self.organism == 'animal':
            logger.info("Running SecretomeP v1.0h to predict non-classical signal peptides.")
            self.seqs_separation()
            multiple_secretome(self.secretomep, self.sub_file_list, self.threads)
            s_dict = defaultdict()
            for files in self.sub_file_list:
                result_file = files[1]
                with open(result_file, 'r') as f:
                    for each_line in f.readlines()[5:]:
                        if '#' not in each_line:
                            s_list = each_line.strip().split('\t')
                            protein = s_list[0].strip()
                            n_score = s_list[1].strip()
                            s_dict[protein] = n_score
            for each_protein in self.protein_dict.keys():
                if each_protein in s_dict:
                    each_n_score = s_dict[each_protein]
                else:
                    each_n_score = '0.000'
                self.protein_dict[each_protein].append(each_n_score)
        else:
            for each_protein in self.protein_dict.keys():
                each_n_score = '0.000'
                self.protein_dict[each_protein].append(each_n_score)

    def classical_secretion(self):
        """
        此函数用于根据SignalP和Uniprot keywords来自动判断经典分泌蛋白并绘制两种预测方法的Venn图
        """
        logger.info('Predicting classical secreted proteins using SignalP and Uniprot keywords.')
        for protein, info_list in self.protein_dict.items():
            signalp_result = info_list[2]
            if signalp_result != 'OTHER':
                self.classical_dict[protein].append('I')
            keywords = info_list[1]
            if 'Signal' in keywords:
                self.classical_dict[protein].append('II')
        if not self.no_uniprot:
            matplotlib.use('Agg')
            rcParams.update({'figure.autolayout': True})
            common_number = 0
            type_1_number = 0
            type_2_number = 0
            for each_protein, type_list in self.classical_dict.items():
                type_list = list(set(type_list))
                if len(type_list) > 1:
                    common_number += 1
                else:
                    if 'I' in type_list:
                        type_1_number += 1
                    else:
                        type_2_number += 1
            venn2(subsets=(type_1_number, type_2_number, common_number),
                  set_labels=('SignalP', 'Uniprot Keywords'))
            plt.title("Prediction of classical secreted proteins")
            image_pdf = os.path.join(output_dir, 'Classical_secreted_proteins_Venn_diagram.pdf')
            image_png = os.path.join(output_dir, 'Classical_secreted_proteins_Venn_diagram.png')
            plt.savefig(image_pdf)
            plt.savefig(image_png, dpi=500)
            plt.close()

    def non_classical_secretion(self):
        """
        此函数用于根据GO CC和SecretomeP中的NN-score来判断非经典分泌蛋白
        """
        logger.info('Predicting remaining non-classically secreted proteins using SecretomeP and GO CC.')
        for protein, info_list in self.protein_dict.items():
            if protein not in self.classical_dict:
                go_cc = info_list[0]
                if 'extracellular' in go_cc:
                    self.non_classical_dict[protein].append('I')
                elif 'intracellular' not in go_cc:
                    if self.organism == 'animal':
                        nn_score = float(info_list[3])
                        if nn_score > 0.5:
                            self.non_classical_dict[protein].append('II')

    def combine_secretion(self):
        """
        此函数用于整合经典和非经典分泌蛋白的结果
        """
        logger.info('Combining both classical and non-classical secreted proteins.')
        self.temp_result_file = os.path.join(self.temp_dir, 'temp_results.txt')
        result_lines = 'Protein accession\tProtein description\tGene name\tSecreted type\t' \
                       'GO CC\tUniprot Keywords\tSignalP prediction\tNN_score (SecretomeP)\n'
        for protein, info_list in self.protein_dict.items():
            result_list = self.annotation_dict[protein]
            if protein in self.classical_dict:
                result_list.append('Classical')
                result_list.extend(info_list)
                result_lines += '\t'.join(result_list) + '\n'
            elif protein in self.non_classical_dict:
                result_list.append('Non-classical')
                result_list.extend(info_list)
                result_lines += '\t'.join(result_list) + '\n'
            else:
                result_list.append('Other')
                result_list.extend(info_list)
                result_lines += '\t'.join(result_list) + '\n'
        with open(self.temp_result_file, 'w') as o:
            o.write(result_lines)
        classical_number = len(self.classical_dict)
        non_classical_number = len(self.non_classical_dict)
        self.statistics_dict['Classical'] = [classical_number]
        self.statistics_dict['Non-classical'] = [non_classical_number]
        logger.info('Classical and non-classical protein number: [{}] & [{}].'.format(classical_number,
                                                                                      non_classical_number))

    def write_excel(self):
        """
        此函数用于将分泌蛋白预测结果输出为Excel文件并绘制统计条形图
        """
        logger.info('Writing all prediction results to an Excel format file.')
        stat_df = pd.DataFrame(data=self.statistics_dict)
        matplotlib.use('Agg')
        rcParams.update({'figure.autolayout': True})
        ax = stat_df.T.plot.bar(rot=0, legend=None, figsize=(3.6, 4), width=0.3)
        for pl in ax.patches:
            ax.annotate(pl.get_height(), (pl.get_x() + pl.get_width() / 2, pl.get_height() * 1.001),
                        ha='center', va='bottom')
        # for pl in ax.patches:
        #     ax.annotate(str(pl.get_height()), (pl.get_x() * 1.005, pl.get_height() * 1.005))
        plt.title("Prediction of secreted proteins")
        image_pdf = os.path.join(output_dir, 'Secreted_proteins_prediction.pdf')
        image_png = os.path.join(output_dir, 'Secreted_proteins_prediction.png')
        plt.savefig(image_pdf)
        plt.savefig(image_png, dpi=500)
        plt.close()
        df = pd.read_csv(self.temp_result_file, sep='\t', index_col=False)
        if self.no_uniprot:
            df = df.drop('Uniprot Keywords', axis=1)
        if self.organism != 'animal':
            df = df.drop('NN_score (SecretomeP)', axis=1)
        new_df = df.style.set_properties(**{
            'background-color': '#D9D9D9',
            'font-size': '10pt',
            'font-family': 'Times New Roman'})
        writer = pd.ExcelWriter(self.excel_file, engine='xlsxwriter')
        workbook = writer.book
        header_format = workbook.add_format({
            'bold': True,
            'text_wrap': True,
            'center_across': True,
            'align': 'center',
            'valign': 'vcenter',
            'font_size': 10,
            'bottom': 2,
            'top': 2,
            'font_name': 'Times New Roman',
            'fg_color': '#00CD00'})
        font_format = workbook.add_format({'font_name': 'Times New Roman', 'font_size': 10})
        stat_df.to_excel(writer, sheet_name='Summary', index=False)
        worksheet = writer.sheets['Summary']
        worksheet.set_column('A:A', 10, font_format)
        worksheet.set_column('B:B', 10, font_format)
        for col_num, value in enumerate(stat_df.columns.values):
            worksheet.write(0, col_num, value, header_format)
        new_df.to_excel(writer, sheet_name='Prediction_result', index=False)
        worksheet = writer.sheets['Prediction_result']
        worksheet.set_column('A:A', 10, font_format)
        worksheet.set_column('B:B', 40, font_format)
        worksheet.set_column('C:C', 10, font_format)
        worksheet.set_column('D:D', 15, font_format)
        worksheet.set_column('E:E', 40, font_format)
        worksheet.set_column('F:F', 40, font_format)
        worksheet.set_column('G:G', 15, font_format)
        worksheet.set_column('H:H', 10, font_format)
        for col_num, value in enumerate(new_df.columns.values):
            worksheet.write(0, col_num, value, header_format)
        writer.save()
        # os.remove(self.temp_result_file)

    def main(self):
        """
        此函数用于自动运行全部分泌蛋白预测流程
        """
        self.load_data()
        self.fetch_protein_sequences()
        self.run_signalp()
        self.classical_secretion()
        self.run_secretomep()
        self.non_classical_secretion()
        self.combine_secretion()
        self.write_excel()
        shutil.rmtree(self.temp_dir)


if __name__ == '__main__':
    # Run as a script
    __version__ = '1.0'
    start_time = time.time()
    my_path = os.getcwd()
    args = parse_cmdline()
    # Set up logging
    logger, formatter = new_logger()
    # Report arguments, if verbose.
    logger.info('Command line: [%s]' % ' '.join(sys.argv))
    logger.info(args)
    output_dir = args.output
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    p = MainPipeline()
    p.main()
    end_time = time.time()
    logger.info("All jobs have been done: %s." % time.asctime())
    logger.info("Total time taken: %.2fs." % (end_time - start_time))
