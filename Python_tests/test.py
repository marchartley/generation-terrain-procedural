import os
import subprocess
import time


server_commands = []
client_commands = ["Hello\n"]*1

# mainP = subprocess.Popen(["/home/marc/codes/IG3/MultiProcessing/multiS", "5000", "tcp"])
# for command in server_commands:
#     stdout_data = mainP.communicate(input=command)[0]



processes = []
for i in range(1):
    p = subprocess.Popen(["/home/marc/codes/IG3/MultiProcessing/multiC", "5000", "tcp"], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    processes.append(p)


for command in client_commands:
    for i in range(len(processes)):
        p = processes[i]
        p.stdin.write(command)
        p.stdin.flush()

print("Ending")