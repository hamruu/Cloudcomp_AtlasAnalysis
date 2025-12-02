#CONTROLLER FILE
import numpy as np # for numerical calculations such as histogramming
import matplotlib.pyplot as plt # for plotting
#matplotlib_inline.backend_inline.set_matplotlib_formats('pdf', 'svg') # to make plots in pdf (vector) format
from matplotlib.ticker import AutoMinorLocator # for minor ticks
import uproot # for reading .root files
import awkward as ak # to represent nested data in columnar format
import vector # for 4-momentum calculations
import time # for printing time stamps
import aiohttp
import requests # for file gathering, if needed
import pika
import json
import os

import atlasopenmagic as atom
atom.available_releases()
atom.set_release('2025e-13tev-beta')

#Definitions for plotting and root tree opening
skim = "exactly4lep"

defs = {
    r'Data':{'dids':['data']},
    r'Background $Z,t\bar{t},t\bar{t}+V,VVV$':{'dids': [410470,410155,410218,
                                                        410219,412043,364243,
                                                        364242,364246,364248,
                                                        700320,700321,700322,
                                                        700323,700324,700325], 'color': "#124405" }, # purple
    r'Background $ZZ^{*}$':     {'dids': [700600],'color': "#2bcb1c" },# red
    r'Signal ($m_H$ = 125 GeV)':  {'dids': [345060, 346228, 346310, 346311, 346312,
                                          346340, 346341, 346342],'color': "#43CCF1" },# light blue
}

samples   = atom.build_dataset(defs, skim=skim, protocol='https', cache=True)

#Setting up messaging with rabbit mq, allowing controller and worker to communicate
rabbit_host = os.environ["RABBIT_HOST"]
params = pika.ConnectionParameters(host=rabbit_host, heartbeat = 600) #Avoids rabbitmq timeout and worker shutdown when processing large data sets or on bad wifi
connection = pika.BlockingConnection(params)
channel = connection.channel()

channel.queue_declare(queue= "taskqueue")
channel.queue_declare(queue = "resultqueue")

all_data = {} #To accumulate the hist data that comes back from the workers

#Sending off root files as tasks for worker nodes
tasks_sent = 0
results_received = 0

for s in samples:

    #Print which sample is being processed
    print('Processing '+s+' samples')

    #Loop over each file
    for val in samples[s]['list']:
        tasks_sent += 1

        task = {"sample": s, "filepath": val}
        channel.basic_publish(exchange='',
            routing_key='taskqueue',
            body=json.dumps(task)
        )

def callback(channel, method, properties, body):
    global results_received
    results_received += 1

    result = json.loads(body.decode())
    sample = result["sample"]
    result = result["result"]
    #Constructing the full data set from worker fragments
    if sample not in all_data:
        all_data[sample] = np.array(result)
    else:
        all_data[sample] += np.array(result)



    channel.basic_ack(delivery_tag=method.delivery_tag)

    if results_received == tasks_sent:
        print("All results received.")
        print("Press Ctrl + C, followed by docker compose down to close service.")
        channel.stop_consuming()
        connection.close()
    

channel.basic_consume(queue="resultqueue", on_message_callback=callback)
channel.start_consuming()

#Plotting code
xmin = 80
xmax = 250
step = 2.5

bin_edges = np.arange(xmin, xmax + step, step)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

data = all_data["Data"]

fig, ax = plt.subplots(figsize=(12, 8))

#Data with errors
ax.errorbar(
    bin_centers,
    data,
    yerr=np.sqrt(data),
    fmt='o',
    color='black',
    label='Data'
)

#MC samples
for s in samples:
    if s == "Data":
        continue
    ax.bar(
        bin_centers,
        all_data[s],
        width=step,
        color=defs[s]["color"],
        alpha=0.6,
        label=s
    )

ax.set_xlabel(r'4-lepton invariant mass $\mathrm{m_{4l}}$ [GeV]',
                    fontsize=13, x=1, horizontalalignment='right' )
ax.set_ylabel('Events / '+str(step)+' GeV',
                        y=1, horizontalalignment='right')

# Add text 'ATLAS Open Data' on plot
plt.text(0.2, # x
            0.93, # y
            'ATLAS Open Data', # text
            transform=ax.transAxes, 
            fontsize=16 )

# Add text 'for education' on plot
plt.text(0.2, # x
            0.88, # y
            'for education', # text
            transform=ax.transAxes, 
            style='italic',
            fontsize=12 )

# Add energy and luminosity
lumi_used = str(36.6) # luminosity to write on the plot
plt.text(0.2, # x
            0.82, # y
            r'$\sqrt{s}$=13 TeV,$\int$L dt = '+lumi_used+' fb$^{-1}$', # text
            transform=ax.transAxes,fontsize=16 )

# Add a label for the analysis carried out
plt.text(0.2, # x
            0.76, # y
            r'$H \rightarrow ZZ^* \rightarrow 4\ell$', # text
            transform=ax.transAxes,fontsize=16 )
ax.legend()

plt.tight_layout()
plt.savefig("/plots/greenHZZplot.png", dpi=150)
plt.close()
