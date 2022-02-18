import setuptools
import os


setuptools.setup(

	name="amplicon-seq",

	version="0.0.4",
	packages=["amplicon_seq"],
	license="MIT",
	long_description="Amplicon sequencing command line tool",
	scripts= ["scripts/%s" % x for x in os.listdir("scripts")],
)
