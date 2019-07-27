# codevp
Ray Tracing and Optical Design in Python

CodeVP is a sequential ray tracing package written in Python for designing optical systems in the geometric optics regime. It features lens optimization from Scipy and (very soon) Scikit-learn. CodeVP is currently in active development and any collaborators would be welcome.

To use CodeVP I suggest cloning the repository and using the command "pip3 install ." in the downloaded directory. You can also download CodeVP directly through

"pip3 install codevp" (outdated most likely)

but the pip version is likely out of date since this software is in active development and features are being frequently added, most noteably lens optimization which is not (yet) available in the pypi version. To get started using the software, I suggest looking at the examples provided in the example_projects, optimization_examples, and validation files. The UI for CodeVP most parallel's BEAM4, with CodeVP's ray tracing algorithm recreated mostly from Spencer and Murty's iconic paper. Some changes were made to Spencer's algorithm however, and the plotting and optimization software was developed independently by myself.

Contact me at gavinniendorf@gmail.com to get involved in CodeVP's development or for questions on how to use it.

- Gavin Niendorf
