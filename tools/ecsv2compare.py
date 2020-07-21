import csv, os, re
import logging as log
import argparse
import xlsxwriter as xlwt


def main ():
	options = _set_options()
	log.debug("Merge csv file ")
	wb = xlwt.Workbook(options.out)
	blast_query_id_list = list()
	for filename in options.blast_files:
		(f_path, f_name) = os.path.split(filename)
		(f_short_name_tmp, f_extension) = os.path.splitext(f_name)
		f_short_name = f_short_name_tmp.replace("scaffold.", "")
		ws = wb.add_worksheet(f_short_name)
		spamReader = csv.reader(open(filename, "r"), delimiter="\t")
		for rowx, row in enumerate(spamReader):
			cell_format = wb.add_format({})
			if len(row) > 13 :
				taxonomy = row[14]
				if taxonomy.startswith( 'Virus' ):
					cell_format = wb.add_format({"bg_color":"#CCFFFF"})
				else:
					cell_format = wb.add_format({})
					blast_query_id_list.append(row[1])
			ws.write_row(rowx,0, row, cell_format)
		ws.autofilter('A1:P1')
	if (options.rps_file != ""):
		for filename in options.rps_file:
			pattern = r"(.)*Virus(.)*;"
			regex = re.compile(pattern)
			(f_path, f_name) = os.path.split(filename)
			(f_short_name_tmp, f_extension) = os.path.splitext(f_name)
			f_short_name = f_short_name_tmp.replace("scaffold.", "")
			ws = wb.add_worksheet(f_short_name)
			spamReader = csv.reader(open(filename, "r"), delimiter="\t")
			for rowx, row in enumerate(spamReader):
				cell_format = wb.add_format({})
				if len(row) > 8 :
					taxonomy = row[9] # col J
					if regex.match(taxonomy) is not None:
						cell_format = wb.add_format({"bg_color":"#CCFFFF"}) # row color as blue if Viruses
						if row[0] in blast_query_id_list:
							cell_format = wb.add_format({"bg_color":"#FF6347"}) # row color as red if pfam_query_id not in blast_query_id
					else:
						cell_format = wb.add_format({})
				ws.write_row(rowx,0, row, cell_format)
			ws.autofilter('A1:M1')
	wb.close()


def _set_options ():
	parser = argparse.ArgumentParser()
	parser.add_argument('-b','--blast',help='A csv Blast file.',action='append',required=True,type=str,dest='blast_files')
	parser.add_argument('-r','--rps',help='A csv pfam file.',action='append',required=False,type=str,dest='rps_file')
	parser.add_argument('-o','--out',help='The output file.',action='store',type=str,default='./',dest="out")
	parser.add_argument('-v','--verbosity',help='Verbose level', action='store',type=int,choices=[1,2,3,4],default=1)
	args = parser.parse_args()
	return args

def _set_log_level(verbosity):
	if verbosity == 1:
		log_format = '%(asctime)s %(levelname)-8s %(message)s'
		log.basicConfig(level=log.INFO,format=log_format)
	elif verbosity == 3:
		log_format = '%(filename)s:%(lineno)s - %(asctime)s %(levelname)-8s %(message)s'
		log.basicConfig(level=log.DEBUG,format=log_format)


if __name__ == "__main__":
	main()