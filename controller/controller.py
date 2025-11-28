#This is the controller.

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
###################################################################################
################################ defining variables ###############################
###################################################################################
skim = "exactly4lep"

defs = {
    r'Data':{'dids':['data']},
    r'Background $Z,t\bar{t},t\bar{t}+V,VVV$':{'dids': [410470,410155,410218,
                                                        410219,412043,364243,
                                                        364242,364246,364248,
                                                        700320,700321,700322,
                                                        700323,700324,700325], 'color': "#6b59d3" }, # purple
    r'Background $ZZ^{*}$':     {'dids': [700600],'color': "#ff0000" },# red
    r'Signal ($m_H$ = 125 GeV)':  {'dids': [345060, 346228, 346310, 346311, 346312,
                                          346340, 346341, 346342],'color': "#00cdff" },# light blue
}

samples   = atom.build_dataset(defs, skim=skim, protocol='https', cache=True)

def calc_mass(lep_pt, lep_eta, lep_phi, lep_e):
    p4 = vector.zip({"pt": lep_pt, "eta": lep_eta, "phi": lep_phi, "E": lep_e})
    invariant_mass = (p4[:, 0] + p4[:, 1] + p4[:, 2] + p4[:, 3]).M # .M calculates the invariant mass
    return invariant_mass

lumi = 36.6
def calc_weight(weight_variables, events):
    total_weight = lumi * 1000 / events["sum_of_weights"]
    for variable in weight_variables:
        total_weight = total_weight * abs(events[variable])
    return total_weight
###############################################################################################################
##################### Variables defined #######################################################################
###############################################################################################################

################### message sending and receiving stuff ################
rabbit_host = os.environ["RABBIT_HOST"]
params = pika.ConnectionParameters(host=rabbit_host)
connection = pika.BlockingConnection(params)
channel = connection.channel()

channel.queue_declare(queue= "taskqueue")
channel.queue_declare(queue = "resultqueue")
################# message sending and receiving stuff ##############


all_data = {} #To accumulate the hist data that comes back from the workers

# Analysis start

tasks_sent = 0
results_received = 0

for s in samples:

    # Print which sample is being processed
    print('Processing '+s+' samples')

    # Define empty list to hold data
    frames = []

    # Loop over each file
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
    filepath = result["filepath"]
    hist_segment = result["result"]

    channel.basic_ack(delivery_tag=method.delivery_tag)

    if results_received == tasks_sent:
        print("All results received.")
        channel.stop_consuming()
        connection.close()
    



    

channel.basic_consume(queue="resultqueue",on_message_callback=callback)
channel.start_consuming()
        
#NEEDS TO RECEIVE AND PLOT!
#Read up on rabbitmq docs, persistent mode for delivery, durable queues.