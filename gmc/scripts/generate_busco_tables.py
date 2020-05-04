import os
import glob
import csv
from gmc.scripts.analyse_busco import read_full_table, read_tx2gene, get_busco_categories
from collections import Counter

class BuscoTableGenerator:
	REVIEW_TABLE_HEADER = ["Busco ID", "Transcript ID", "Busco Status", "Coordinates", "prepare TID (C or D)", "prepare TID coordinates"]

	def import_tx2gene_data(self, tx2gene_data_dir):
		self.tx2gene = dict()
		for f in os.listdir(tx2gene_data_dir):
			f = os.path.join(tx2gene_data_dir, f)
			self.tx2gene.update(read_tx2gene(f))

	def import_txcoords(self, txcoords_prepare, txcoords_pick):
		self.txcoords_prepare = dict(line.strip().split("\t") for line in open(txcoords_prepare))
		self.txcoords_pick = dict(line.strip().split("\t") for line in open(txcoords_pick))

	def parse_rundata(self, busco_run_path):
		for d in os.listdir(busco_run_path):
			if d.endswith("_final") or d == "genome":
				f = glob.glob(os.path.join(busco_run_path, d, d, "run_*", "full*"))[0]
				self.run_tables[d], complete_buscos, missing_buscos, fragmented_buscos = read_full_table(f, is_pick=d.endswith("_final"))
				if d == "proteins_final":
					self.complete_busco_proteins_final = complete_buscos
					self.fragmented_busco_proteins[d] = fragmented_buscos
			else:
				for dd in glob.glob(os.path.join(busco_run_path, d, "*")):
					if os.path.basename(dd) != "input":
						f = glob.glob(os.path.join(dd, "run_*", "full*"))[0]
						self.run_tables["{}/{}".format(d, os.path.basename(dd))], complete_buscos, missing_buscos, fragmented_buscos = read_full_table(f, self.tx2gene)
						if os.path.basename(d).startswith("proteins"):
							self.complete_busco_proteins[dd] = complete_buscos
							if not self.missing_busco_proteins:
								self.missing_busco_proteins.update(missing_buscos)
							else:
								self.missing_busco_proteins.intersection_update(missing_buscos)
							self.fragmented_busco_proteins[dd] = fragmented_buscos
						else:
							self.complete_busco_transcripts[dd] = complete_buscos
							if not self.missing_busco_transcripts:
								self.missing_busco_transcripts.update(missing_buscos)
							else:
								self.missing_busco_transcripts.intersection_update(missing_buscos)

	def __init__(self, tx2gene_data_dir, txcoords_prepare, txcoords_pick, busco_run_path):
		self.import_tx2gene_data(tx2gene_data_dir)
		self.import_txcoords(txcoords_prepare, txcoords_pick)
		
		self.run_tables = dict()
		self.complete_busco_proteins, self.complete_busco_transcripts = dict(), dict()
		self.complete_busco_proteins_final = dict()
		self.missing_busco_proteins, self.missing_busco_transcripts = set(), set()
		self.fragmented_busco_proteins = dict()

		self.parse_rundata(busco_run_path)

	def write_review_table(self, prefix):
		review_proteins = set()
		for protein_set in self.complete_busco_proteins.values():
			review_proteins.update(protein_set)
		review_proteins.difference_update(self.complete_busco_proteins_final)
		with open(prefix + ".review_table", "w") as review_out:
			print(*BuscoTableGenerator.REVIEW_TABLE_HEADER, sep="\t", flush=True, file=review_out)
			for bid in sorted(review_proteins):
				tid = ",".join(self.fragmented_busco_proteins["proteins_final"].get(bid, list()))
				busco_status = "fragmented" if tid else "missing"
				tid_coords = self.txcoords_pick.get(tid, None)
				prepare_tids, prepare_coords = list(), list()
				for protein_set in self.complete_busco_proteins.values():
					prepare_tids.extend(item[0] for item in protein_set.get(bid, list()))
				prepare_coords = [self.txcoords_prepare.get(ptid, None) for ptid in prepare_tids]
				row = [bid, tid, busco_status, tid_coords, ",".join(prepare_tids), ",".join(prepare_coords)]
				print(*row, sep="\t", flush=True, file=review_out)

	def write_raw_data(self, prefix): 
		with open(prefix + ".raw", "w") as raw_out:
			for k, v in self.run_tables.items():
				print(k, v, file=raw_out, sep="\t", flush=True)

	def write_busco_table(self, prefix, max_copy_number):
		with open(prefix, "w") as table_out:
			print("Busco Plots", *self.run_tables.keys(), sep="\t", flush=True, file=table_out)
			for cat in get_busco_categories(max_copy_number=max_copy_number):
				cat_lbl = cat
				if cat.startswith("Complete_"):
					copies = cat.split("_")[1]
					cat_lbl = "Complete (single copy)" if copies == "1" else "Complete ({} copies)".format(copies)
				
				print(cat_lbl, *(v[cat] for v in self.run_tables.values()), sep="\t", flush=True, file=table_out)

			complete_busco_proteins_, complete_busco_transcripts_ = set(), set()
			for set_ in self.complete_busco_proteins.values():
				complete_busco_proteins_.update(set_)
			for set_ in self.complete_busco_transcripts.values():
				complete_busco_transcripts_.update(set_)
			print("# best achievable protein BUSCO count: {}".format(len(complete_busco_proteins_)), flush=True, file=table_out)
			print("# best achievable transcript BUSCO count: {}".format(len(complete_busco_transcripts_)), flush=True, file=table_out)
			print("# lowest achievable missing protein BUSCO count: {}".format(len(self.missing_busco_proteins)), flush=True, file=table_out)
			print("# lowest achievable missing transcript BUSCO count: {}".format(len(self.missing_busco_transcripts)), flush=True, file=table_out)
