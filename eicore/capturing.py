from io import StringIO
import sys

class Capturing(list):
	'''This class is intended to capture the stdout and stderr from a process and display it inline with the host process'''
	def __enter__(self):
		self._stdout = sys.stdout
		self._stderr = sys.stderr
		sys.stdout = self._stringioout = StringIO()
		sys.stderr = self._stringioerr = StringIO()
		return self

	def __exit__(self, *args):
		self.extend(self._stringioout.getvalue().splitlines())
		self.extend(self._stringioerr.getvalue().splitlines())
		del self._stringioout  # free up some memory
		del self._stringioerr  # free up some memory
		sys.stdout = self._stdout
		sys.stderr = self._stderr