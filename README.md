# nursery_shelling_labels
Python code to generate plant breeding nursery inventories and shelling/threshing labels for seed packets.

A script is used to combine input file information, add additionl information, then call the createLabels() function of the nursery module:

![](https://github.com/ncsumaize/nursery_shelling_labels/blob/master/images/Workflow_main_shelling_label_script.png)

The createLabels() function is in the nursery module and it calls a bunch of other functions within that module to produce the new harvest labels:

![](https://github.com/ncsumaize/nursery_shelling_labels/blob/master/images/Workflow_nursery_module.png)
  
[Sample nursery info file used for input](nursery_sample_info.csv) 

[Sample nursery harvest file used for input](harvest_notes_sample.csv)  

[iPython notebook to combine the info and harvest information and call the nursery functions](Nursery_sample.ipynb)  
  
[nursery module. Save as nursery.py in your working directory](nursery.py)
