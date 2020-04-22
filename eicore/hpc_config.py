import json


class HpcConfig:
	def __init__(self, cfg):
		self.__records = json.load(open(cfg))
		self.default_cfg = self.__records["__default__"]


	def __get_resource(self, rulename, resource):
		return self.__records.get(rulename, dict()).get(resource, self.default_cfg[resource])

	def get_cores(self, rulename):
		return int(self.__get_resource(rulename, "cores"))

	def get_memory(self, rulename):
		return int(self.__get_resource(rulename, "memory"))

	def get_time(self, rulename):
		return self.__get_resource(rulename, "time")

	def get_exclude_list(self, rulename):
		exclude = self.__get_resource(rulename, "exclude")
		return "-x {}".format(exclude) if exclude else ""
