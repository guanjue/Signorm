### get parameters
args = commandArgs(trailingOnly=TRUE)

input_file_list = args[1]

input_file_names = read.table(input_file_list, header = FALSE)

for (files in input_file_names){
	print(files[1])
	print(files[2])
}

