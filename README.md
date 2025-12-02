## **Cloud Computing Project – SCIFM0004**

For this project, a system was built to analyse atlas experiment data from the large hadron collider.

The system consists of three components:

### **1. Controller**
Sends file paths to the worker nodes. Then, gathers back data from the workers and combines it all into one plot.

### **2. Worker(s)**
Carries out the actual processing on experimental and monte carlo data.

### **3. RabbitMQ**
Allows the controller and worker(s) to communicate with eachother.

---

## **How to Run**


```
docker compose up --scale worker=N
```

Where N is the desired number of workers.


This project is based off a pre-existing atlas jupyter notebook:  https://github.com/atlas-outreach-data-tools/notebooks-collection-opendata/blob/master/13-TeV-examples/uproot_python/HZZAnalysis.ipynb
