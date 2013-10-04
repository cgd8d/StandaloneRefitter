
BaseDir = "/scratch1/scratchdirs/claytond/ChargeAndLight"
NumEntriesPerJob = 10
NumJobs = 4

for i in range(NumJobs):
    with open(BaseDir + "/InOutFiles/infile%04i.txt" % i, 'w') as infile:
        infile.write("proc00004544-000.root\nrun00004544-000.root\n" + BaseDir + "/DenoisedFiles/Denoised%04i.root\n" % i + str(i*NumEntriesPerJob) + "\n" + str(NumEntriesPerJob) + "\n1\n")

