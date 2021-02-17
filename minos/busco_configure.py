import os
import glob

BUSCO_LEVELS = {
	"proteins", "proteome", "prot", "p",
	"transcripts", "transcriptome", "tran", "t",
	"genome", "geno", "g",
	"all", "none", "off"
}


class BuscoConfiguration(dict):
	@staticmethod
	def parse_busco_levels(levels):
		if not levels:
			return False, False, False  #  while developing
			# return True, False, False
		levels = set(l.lower() for l in levels.split(","))
		invalid_levels = levels.difference(BUSCO_LEVELS)
		if invalid_levels:
			raise ValueError(
				"Invalid busco levels specified with --busco-level option ({}). Valid levels are {}.".format(invalid_levels, BUSCO_LEVELS))
		switch_off = {"off", "none"}.intersection(levels)
		unique_levels = set(level[0] for level in levels.difference({"off", "none"}))
		do_proteins = unique_levels.intersection({"a", "p"}) and not switch_off
		do_transcripts = unique_levels.intersection({"a", "t"}) and not switch_off
		do_genome = unique_levels.intersection({"a", "g"}) and not switch_off
		return tuple(map(bool, (do_proteins, do_transcripts, do_genome)))

	def _check_precomputed_genome_run(self, bg_run):
		precomputed_genome = None
		busco_data = {
			"summary": {"fn": "*short_summary*.txt", "label": "genome summary"},
			"full_table": {"fn": "full_table.tsv", "label": "genome full table"},
			"missing_busco_list": {"fn": "missing_busco_list.tsv", "label": "missing busco list"}
		}
		if bg_run is not None:
			precomputed_genome = dict()
			for key, vals in busco_data.items():
				try:
					precomputed_genome[key] = glob.glob( os.path.join(bg_run, vals["fn"]))[0]
				except:
					raise ValueError("No valid {label} file in {bg_run}".format(label=vals["label"], bg_run=bg_run))
		return precomputed_genome

	def __init__(self, args):
		self["proteins"], self["transcriptome"], self["genome"] = BuscoConfiguration.parse_busco_levels(args.busco_level)
		self["proteins"] = self["proteins"] or args.busco_scoring is not None
		if self["genome"] and args.busco_genome_run:
			raise ValueError("BUSCO genome was selected ({}) together with --busco_genome_run ({}).".format(args.busco_level, args.busco_genome_run))
		if any(self.values()) and (args.busco_lineage is None or not os.path.exists(args.busco_lineage)):
			raise ValueError("BUSCO analysis requested (P:{}, T:{}, G:{}) but no valid lineage specified ({}).".format(
				*self.values(),
				args.busco_lineage
			))
		if args.busco_lineage and not set(args.busco_level.split(",")).issubset(BUSCO_LEVELS.difference({"off", "none"})):
			raise ValueError("BUSCO lineage (--busco-lineage ({})) requires --busco-level ({} - user provided) not in {{off,none}})".format(
				args.busco_lineage, args.busco_level))

		self["precomputed_genome"] = self._check_precomputed_genome_run(args.busco_genome_run)
		self["lineage"] = args.busco_lineage

	def to_dict(self):
		return {k: v for k, v in self.items()}
