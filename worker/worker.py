import pika
rabbit_host = os.environ["RABBIT_HOST"]
params = pika.ConnectionParameters(host=rabbit_host)
connection = pika.BlockingConnection(params)
channel = connection.channel()

channel.queue_declare(queue= "task_queue")

messages_received = 0
#DEFINE PROCESS FILE AND PUT IT IN CALLBACK!!

def process_file(val, s): #Val and S have come from channel queue consume and so can be used here.
     # start the clock
    start = time.time()
    print("\t"+val+":")
    fileString = val

        # Open file
    tree = uproot.open(fileString + ":analysis")

    sample_data = []


    for data in tree.iterate(variables + weight_variables + ["sum_of_weights", "lep_n"], library="ak", entry_stop=tree.num_entries*fraction):
        nIn = len(data)

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

            # Number Cuts

            # Lepton cuts

        lep_type = data['lep_type']
        data = data[~cut_lep_type(lep_type)]
        lep_charge = data['lep_charge']
        data = data[~cut_lep_charge(lep_charge)]

            # Invariant Mass
        data['mass'] = calc_mass(data['lep_pt'], data['lep_eta'], data['lep_phi'], data['lep_e'])

            # Store Monte Carlo weights in the data
        if 'data' not in s: # Only calculates weights if the data is MC
            data['totalWeight'] = calc_weight(weight_variables, data)

            # Append data to the whole sample data list
        sample_data.append(data)

        if not 'data' in val:
            nOut = sum(data['totalWeight']) # sum of weights passing cuts in this batch
        else:
            nOut = len(data)

        elapsed = time.time() - start # time taken to process
        print("\t\t nIn: "+str(nIn)+",\t nOut: \t"+str(nOut)+"\t in "+str(round(elapsed,1))+"s") # events before and after




def callback(channel, method, properties, body):
      messages_received+=1
      print (f"Received message {messages_received}")

      #PROCESS FILE HERE!!!

channel.basic_consume(queue="task_queue",
                      on_message_callback=callback)


