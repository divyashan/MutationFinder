from pdf2txt import main
import os

gene_name = 'CFTR'

def process_pdf_files(gene_name, pdf_path):
	home = os.path.join("/Users", "divya")
	papers_by_gene = os.path.join( home, "Google Drive", "Curation Papers", "Curation Papers by Gene", gene_name)
	if not os.path.exists(gene_name+"_text_files"):
		os.makedirs(gene_name + '_text_files')
	for paper in os.listdir(papers_by_gene):
		with open(os.path.join(papers_by_gene, paper), 'r') as pdfFile:
			print paper
			try:
				pmid = int(paper.split("_")[0])
			except ValueError:
				continue

			output_path = paper[:-4] + '.txt'
			file_path = os.path.join(papers_by_gene, paper)
			try:
				convert(['./pdf2txt.py','-o', output_path, file_path])
			except:
				print "Text extraction not allowed :("

