#WORKER FILE
import numpy as np # for numerical calculations such as histogramming
import uproot # for reading .root files
import awkward as ak # to represent nested data in columnar format
import vector # for 4-momentum calculations
import pika
import json
import os
import pika

#Setting up messaging with rabbit mq, allowing controller and worker to communicate
rabbit_host = os.environ["RABBIT_HOST"]
params = pika.ConnectionParameters(host=rabbit_host, heartbeat = 600)
connection = pika.BlockingConnection(params)
channel = connection.channel()

channel.queue_declare(queue= "taskqueue")
channel.queue_declare(queue = "resultqueue")
messages_received = 0

#Defining variables and functions needed for process file
weight_variables = ["filteff","kfac","xsec","mcWeight","ScaleFactor_PILEUP", "ScaleFactor_ELE", "ScaleFactor_MUON", "ScaleFactor_LepTRIGGER"]
variables = ['lep_pt','lep_eta','lep_phi','lep_e','lep_charge','lep_type','trigE','trigM','lep_isTrigMatched',
            'lep_isLooseID','lep_isMediumID','lep_isLooseIso','lep_type']

def cut_trig_match(lep_trigmatch):
    trigmatch = lep_trigmatch
    cut1 = ak.sum(trigmatch, axis=1) >= 1
    return cut1

def cut_trig(trigE,trigM):
    return trigE | trigM

def ID_iso_cut(IDel,IDmu,isoel,isomu,pid):
    thispid = pid
    return (ak.sum(((thispid == 13) & IDmu & isomu) | ((thispid == 11) & IDel & isoel), axis=1) == 4)

def cut_lep_type(lep_type):
    sum_lep_type = lep_type[:, 0] + lep_type[:, 1] + lep_type[:, 2] + lep_type[:, 3]
    lep_type_cut_bool = (sum_lep_type != 44) & (sum_lep_type != 48) & (sum_lep_type != 52)
    return lep_type_cut_bool # True means we should remove this entry (lepton type does not match)

def cut_lep_charge(lep_charge):
    # first lepton in each event is [:, 0], 2nd lepton is [:, 1] etc
    sum_lep_charge = lep_charge[:, 0] + lep_charge[:, 1] + lep_charge[:, 2] + lep_charge[:, 3] != 0
    return sum_lep_charge # True means we should remove this entry (sum of lepton charges is not equal to 0)

def calc_weight(weight_variables, events):
    lumi = 36.6
    total_weight = lumi * 1000 / events["sum_of_weights"]
    for variable in weight_variables:
        total_weight = total_weight * abs(events[variable])
    return total_weight

def calc_mass(lep_pt, lep_eta, lep_phi, lep_e):
    p4 = vector.zip({"pt": lep_pt, "eta": lep_eta, "phi": lep_phi, "E": lep_e})
    invariant_mass = (p4[:, 0] + p4[:, 1] + p4[:, 2] + p4[:, 3]).M # .M calculates the invariant mass
    return invariant_mass


#Process file function that performs the analysis
def process_file(s, val): #Val and S have come from channel queue consume and so can be used here.
    
    fraction = 1

    print("\t"+val+":")
    fileString = val

    xmin = 80
    xmax = 250
    step_size = 2.5
    bin_edges = np.arange(start=xmin, # The interval includes this value
                        stop=xmax+step_size, # The interval doesn't include this value
                        step=step_size ) # Spacing between values
    
    #Open file
    tree = uproot.open(fileString + ":analysis")

    frames = np.zeros(len(bin_edges)-1)
    

    for data in tree.iterate(variables + weight_variables + ["sum_of_weights", "lep_n"], library="ak", entry_stop=tree.num_entries*fraction):

        data = data[cut_trig(data.trigE, data.trigM)]
        data = data[cut_trig_match(data.lep_isTrigMatched)]

            # Record transverse momenta (see bonus activity for explanation)
        data['leading_lep_pt'] = data['lep_pt'][:,0]
        data['sub_leading_lep_pt'] = data['lep_pt'][:,1]
        data['third_leading_lep_pt'] = data['lep_pt'][:,2]
        data['last_lep_pt'] = data['lep_pt'][:,3]

            # Cuts on transverse momentum
        data = data[data['leading_lep_pt'] > 20]
        data = data[data['sub_leading_lep_pt'] > 15]
        data = data[data['third_leading_lep_pt'] > 10]

        data = data[ID_iso_cut(data.lep_isLooseID,
                                   data.lep_isMediumID,
                                   data.lep_isLooseIso,
                                   data.lep_isLooseIso,
                                   data.lep_type)]

        #Cuts
        lep_type = data['lep_type']
        data = data[~cut_lep_type(lep_type)]
        lep_charge = data['lep_charge']
        data = data[~cut_lep_charge(lep_charge)]

        #Calculation of invariant mass
            # Invariant Mass
        data['mass'] = calc_mass(data['lep_pt'], data['lep_eta'], data['lep_phi'], data['lep_e'])

            # Store Monte Carlo weights in the data
        if 'data' not in s.lower(): # Only calculates weights if the data is MC
            data['totalWeight'] = calc_weight(weight_variables, data)
        else:
            data['totalWeight'] = 1

        #Make hist for this loop
        loopwise_hist, _ = np.histogram(ak.to_numpy(data['mass']), bins=bin_edges, weights= ak.to_numpy(data['totalWeight']) )
        #Accumulate the data to frames
        frames += loopwise_hist
    #Return the data, send hists back.
    return frames.tolist()    


 




def callback(channel, method, properties, body):
    global messages_received
    messages_received+=1
    print (f"Received message {messages_received}")
    
    task = json.loads(body.decode())
    sample = task["sample"]
    filepath = task["filepath"]

    data = process_file(sample, filepath) #Processes data for filepath as it is received
    
    result = {"sample": sample, "filepath": filepath , "result": data}
    channel.basic_publish(exchange='',
                          routing_key='resultqueue',
                          body=json.dumps(result)) #Results sent back to controller node
    
    channel.basic_ack(delivery_tag=method.delivery_tag)


channel.basic_consume(queue="taskqueue",
                      on_message_callback=callback)

channel.start_consuming()
