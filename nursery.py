# -*- coding: cp1252 -*-
from pandas import Series, DataFrame #we are going to use these a lot
import pandas as pd
import os as os
import re as re
import math as math
from time import strftime 

#####################################################################################
#DEFINE ALL NEEDED FUNCTIONS
#####################################################################################
#Need to create a column called 'Row_ID' that has full information on the row including year and location.
#Eg, row 1 in 14 CL should have Row_ID = 14CL0001. We use the nusery_pre string identified at beginning to define this for each nusery
def rowToRowID(r):
    if r['Row'] < 10: 
        Row_ID = nursery_pre + '000' + str(r['Row'])
    elif r['Row'] >= 10 and r['Row'] < 100:
        Row_ID  = nursery_pre + '00' + str(r['Row'])
    elif r['Row'] >= 100 and r['Row'] < 1000:
        Row_ID  = nursery_pre + '0' + str(r['Row'])
    else: Row_ID  = nursery_pre + str(r['Row'])
    return Row_ID


#Here are the patterns we need to deal with in going from plot generation to seed generation:

#Plot generation:     
#F1
#Fx
#Fx:y
#S0
#Sx
#Sx:y
#BCzF1
#BCzFx:y
#Inbred
#Inbred?
#RMx  (or RMwSx:ySbz)
 #
#Fx:ySibz
 
#Poll_Type:
#SF: selfing
#SB: sibbing within the row
#CR: crossing between paired rows
#TC: topcross for paired rows
#BC: backcross for paired rows
#SF/CR: selfing some plants, crossing others 
#SF/TC: selfing some plants, topcrossing others 
#SF/BC: selfing some plants, backcrossing others 
#FC: free-cross among arbitrary set of rows
#IN: intermate among arbitrary set of rows
#Function to generate seed generations from selfing, note that sibbing generation does not matter for this:

def selfing(r, bc_gen, rm_gen, prefix_char, prefix_num, suffix, earno_self, earq_self):
    #Bulk shelling
    if r['shell'] == 'BULK': 
        if bc_gen:           
            if suffix: seed_gen =  'BC' + str(bc_gen) + prefix_char + str(prefix_num) + ':' + str(suffix + 1) #if suffix is non-empty
            else: seed_gen =  'BC' + str(bc_gen) + prefix_char + str(prefix_num + 1) #if suffix is ''
        elif rm_gen:
            if suffix: seed_gen =  'RM' + str(rm_gen) + prefix_char + str(prefix_num) + ':' + str(suffix + 1) #if suffix is non-empty
            else: seed_gen =  'RM' + str(rm_gen) + prefix_char + str(prefix_num + 1) #if suffix is ''         
        else:
            if suffix: seed_gen =  prefix_char + str(prefix_num) + ':' + str(suffix + 1) #if suffix is non-empty
            elif prefix_char == 'F' and prefix_num == 1: #special case of selfing from an F1
                seed_gen = 'F2'
            else: seed_gen =  prefix_char + str(prefix_num + 1) #if suffix is ''
 
        #some lines are considered 'Inbred' because they were obtained as inbreds
        if prefix_char == 'Inbred': seed_gen =  'Inbred'
        if prefix_char == 'Inbred?': seed_gen =  'Inbred?'
        
                
    #Generate seed_ID for case of bulk shelling
        seed_ID = r['Row_ID'] + '-BLK'
        
        return DataFrame({'seed_ID' : seed_ID,
                          'Material' : r['Material'],
                          'seed_gen' : seed_gen,
                          'ear_number': earno_self,
                          'ear_quality' : earq_self,
                          'source_ID' : r['Source ID']}, index = [0])
 
    
    #Single ear shelling        
    else: 
        if bc_gen:
            if suffix: seed_gen =  'BC' + str(bc_gen) + prefix_char + str(suffix) + ':' + str(suffix + 1) #if suffix is non-empty
            else: seed_gen =  'BC' + str(bc_gen) + prefix_char + str(prefix_num) + ':' + str(prefix_num + 1) #if suffix is ''         
        elif rm_gen:
            if suffix: seed_gen =  'RM' + str(rm_gen) + prefix_char + str(suffix) + ':' + str(suffix + 1)  #if suffix is non-empty
            else: seed_gen =  'RM' + str(rm_gen) + str(prefix_num) + ':' + str(prefix_num + 1) #if suffix is ''         
        else:
            if suffix: seed_gen =  prefix_char + str(suffix) + ':' + str(suffix + 1) #if suffix is non-empty
            elif prefix_char == 'F' and prefix_num == 1: #special case of selfing from an F1
                seed_gen = 'F2'
            else: seed_gen =  prefix_char + str(prefix_num) + ':' + str(prefix_num + 1) #if suffix is ''   
 
        #some lines are considered 'Inbred' because they were obtained as inbreds
        if prefix_char == 'Inbred': seed_gen =  'Inbred'
        if prefix_char == 'Inbred?': seed_gen =  'Inbred?'
                
        result = DataFrame()
        for ear in range(1, int(earno_self) + 1):
            seed_ID = r['Row_ID'] + '-' + str(ear)
            ear_to_row = DataFrame({'seed_ID' : seed_ID,
                                    'Material' : r['Material'],
                                    'seed_gen' : seed_gen,
                                    'ear_quality' : earq_self,
                                    'source_ID' : r['Source ID']}, index = [0])
            result = result.append(ear_to_row)
        return result
 
 
#Function to generate seed generations from sibbing, note that we assume sibbing is followed by bulk shelling as a method of maintaining seed generation, rather than inbreeding:
def sibbing(r, bc_gen, rm_gen, prefix_char, prefix_num, suffix, sib, earno_self, earq_self):
    
    if r['shell'] != 'BULK':  raise Exception("Sibbing without bulking not currently allowed! See function sibbing()")
        
    if bc_gen:
        if suffix: 
            if sib: seed_gen = 'BC' + str(bc_gen) + prefix_char + str(prefix_num) + ':' + str(suffix) + 'SB' + str(sib + 1) #if suffix and sib both non-empty
            else: seed_gen = 'BC' + str(bc_gen) + prefix_char + str(prefix_num) + ':' + str(suffix) + 'SB1'  #if suffix non-empty and sib empty
        else: 
            if sib: seed_gen = 'BC' + str(bc_gen) + prefix_char + str(prefix_num) + 'SB' + str(sib + 1)  #if suffix is '' and sib is non-empty
            else: seed_gen = 'BC' + str(bc_gen) + prefix_char + str(prefix_num) + 'SB1' #if both suffix and sib empty
    
    elif rm_gen:
        if suffix: 
            if sib: seed_gen = 'RM' + str(rm_gen) + prefix_char + str(prefix_num) + ':' + str(suffix) + 'SB' + str(sib + 1) #if suffix and sib both non-empty
            else: seed_gen = 'RM' + str(rm_gen) + prefix_char + str(prefix_num) + ':' + str(suffix) + 'SB1'  #if suffix non-empty and sib empty
        else: 
            if sib: seed_gen = 'RM' + str(rm_gen) + prefix_char + str(prefix_num) + 'SB' + str(sib + 1)  #if suffix is '' and sib is non-empty
            else: seed_gen = 'RM' + str(rm_gen) + prefix_char + str(prefix_num) + 'SB1' #if both suffix and sib empty
            
    else:
        if suffix: 
            if sib: seed_gen = prefix_char + str(prefix_num) + ':' + str(suffix) + 'SB' + str(sib + 1) #if suffix and sib both non-empty
            else: seed_gen = prefix_char + str(prefix_num) + ':' + str(suffix) + 'SB1'  #if suffix non-empty and sib empty
        else: 
            if sib: seed_gen = prefix_char + str(prefix_num) + 'SB' + str(sib + 1)  #if suffix is '' and sib is non-empty
            else: seed_gen = prefix_char + str(prefix_num) + 'SB1' #if both suffix and sib empty
            
    if prefix_char == 'Inbred': seed_gen = 'Inbred'
    if prefix_char == 'Inbred?': seed_gen = 'Inbred?'
                
    seed_ID = r['Row_ID'] + '-BLK'
        
    return DataFrame({'seed_ID' : seed_ID,
                      'Material' : r['Material'],
                      'seed_gen' : seed_gen,
                      'ear_number': earno_self,
                      'ear_quality' : earq_self,
                      'source_ID' : r['Source ID']}, index = [0])
                     
#Helper function to check if a material name includes the cross symbol, and if so if it is surrounded by parentheses
#if it has the cross but no parentheses, it returns the material name wrapped in parentheses
#otherwise, returns the original material string
def parenth(material):
    if type(material) != str:
        print("Material name " + material + " not proper string")
        return(material)
    if re.search(u'\xd7', material):
        if not material.startswith("(") and not material.endswith(")"):
            return("(" + material + ")")
    return(material) #if not already returned

#Function to generate labels from crossing in paired rows
def crossing(r, inputdf, seed_gen, earno_cross, earq_cross):
    #print('inside crossing for row ' + str(r['Row']))

    row = r['Row'] #get the row number of current row
    paired_code = pairs[row] #get the paired_code of current row
    #use dict.get(key, default) to get the paired code of the adjacent rows
    #note that dict[key] returns an exception if key is not in the dict
    #so we have to handle cases where the row is at beginning or end of a series 
    #in such cases, one of the adjacent rows will not exist, and we want to return a paired_code for that row
    #which will not equal the current paired_code, so we just add one to the current paired code as the default for missing key values
    if pairs.get(row - 1, paired_code + 1) == paired_code: pair_row = row - 1 #check if previous row has matching paired_code
    elif pairs.get(row + 1, paired_code + 1) == paired_code: pair_row = row + 1 #else check if next row has matching paired_code
    pair = inputdf[inputdf['Row'] == pair_row] # pull the paired row as a Series from the inputdf data frame

    pair.index = [0] #set the index to zero so we can easily access the contents
    
#    print('inside of crossing function, here are the Source ID values for row and pair:')
#    print(r['Source ID'])
#    print(pair['Source ID'][0])
#    print('inside of crossing function, here are the complete contents of row and pair:')
#    print(r)
#    print(pair)
    
    #use u'\xd7' which is ascii code for the cross symbol, otherwise have encoding problems
    #if either of the parental materials already has a cross symbol but not parenth, it needs to be wrapped in parenthesis
    
    Material = parenth(r['Material']) + u'\xd7' + parenth(pair['Material'][0]) #access the 0 index to get just the Material content itself instead of the Series
    seed_ID = r['Row_ID'] + u'\xd7' + pair['Row_ID'][0][4:]
    source_ID = '(' + r['Source ID'] + ')' + u'\xd7' + '(' + pair['Source ID'][0] + ')'
    
    #Bulk shelling
    if r['shell'] == 'BULK': 
        seed_ID = seed_ID + '-BLK'
        return DataFrame({'seed_ID' : seed_ID,
                          'Material' : Material,
                          'seed_gen' : seed_gen,
                          'ear_number': earno_cross,
                          'ear_quality' : earq_cross,
                          'source_ID' : source_ID}, index = [0])

    else:
        result = DataFrame()
        for ear in range(1, int(earno_cross) + 1):
            ear_ID = seed_ID + '-' + str(ear)
            ear_to_row = DataFrame({'seed_ID' : ear_ID,
                          'Material' : Material,
                          'seed_gen' : seed_gen,
                          'ear_quality' : earq_cross,
                          'source_ID' : source_ID}, index = [0])
            result = result.append(ear_to_row)
        
        return result

#Function to generate backcross labels
def backCrossing(r, inputdf, bc_gen, seed_gen, earno_cross, earq_cross):
    '''
    Need to consider these possibilities:
    1. Making first backcross between inbred and F1 or F1 and inbred
    2. Making later backcross between inbred and BCx or BCx and inbred
    
    '''
    #print('inside backCrossing()')
    #print(r)
    #print('bc_gen = ' + str(bc_gen))
    
    
    row = r['Row'] #get the row number of current row
    paired_code = pairs[row] #get the paired_code of current row
    #use dict.get(key, default) to get the paired code of the adjacent rows
    #note that dict[key] returns an exception if key is not in the dict
    #so we have to handle cases where the row is at beginning or end of a series 
    #in such cases, one of the adjacent rows will not exist, and we want to return a paired_code for that row
    #which will not equal the current paired_code, so we just add one to the current paired code as the default for missing key values
    if pairs.get(row - 1, paired_code + 1) == paired_code: pair_row = row - 1 #check if previous row has matching paired_code
    elif pairs.get(row + 1, paired_code + 1) == paired_code: pair_row = row + 1 #else check if next row has matchin paired_code
    pair = inputdf[inputdf['Row'] == pair_row] # pull the paired row as a Series from the inputdf data frame
    pair.index = [0] #set the index to zero so we can easily access the contents

    #print('pair is')
    #print(pair)
    #print("pair['Gen'][0]")
    #print(str(pair['Gen'][0]))

    Material = parenth(r['Material']) + u'\xd7' + parenth(pair['Material'][0]) #access the 0 index to get just the Material content itself instead of the Series
    #seed_ID = r['Row_ID'] + u'\xd7' + pair['Row_ID'][0][4:] #old version, causes exception
    seed_ID = r['Row_ID'] + u'\xd7' + str(pair['Row_ID'][0][4:])
    source_ID = '(' + r['Source ID'] + ')' + u'\xd7' + '(' + pair['Source ID'] + ')'
    
    #check if paired row is a BC generation 
    match_bc = re.match('BC[0-9]', pair['Gen'][0]) #returns match object only if paired row 'Gen' starts with 'BC'
    #if there was a BC generation designator, we now get the BC generation number, then strip off the BC  prefix and deal with the rest
    if match_bc:
        pair_bc_gen = pair['Gen'][0][2:match_bc.end()] #this gets the bc generation number of the paired row
    else: pair_bc_gen = 0
    #print('Jim added this for debugging!! pair_bc_gen is:')
    #print(pair_bc_gen)
 
    #Case of first backcross of inbred female by F1 male 
    #seed_gen passed from seed_generation() is not correct in this case, need to update it
    if bc_gen == 0 and pair['Gen'][0] == "F1":
        print('backCrossing case 1')
        seed_gen = "BC1F1"
        splitCross = pair['Material'][0].split(u'\xd7') #split the Material into two pieces at the cross symbol
        if r['Material'] == splitCross[0]:
            Material = '(' + r['Material'] + '*2)' + splitCross[1] #eg, (B73*2)Oh43
        elif r['Material'] == splitCross[1]:
            Material = '(' + r['Material'] + '*2)' + splitCross[0] 
        
    #Case of first backcross of F1 female by Inbred male 
    elif bc_gen == 0 and r['Gen'] == "F1":
        print('backCrossing case 2')
        splitCross = r['Material'].split(u'\xd7') #split the Material into two pieces at the cross symbol
        if pair['Material'][0] == splitCross[0]:
            Material = '(' + pair['Material'][0] + '*2)' + splitCross[1] #eg, (B73*2)Oh43
        elif pair['Material'][0] == splitCross[1]:
            Material = '(' + pair['Material'][0] + '*2)' + splitCross[0] 
            
     #Case where first backcross but neither parent has cross symbol, throw an exception
    if (bc_gen == 0) and (not match_bc) and not re.search(u'\xd7', r['Material'] + pair['Material'][0]): 
        print('Inside backCrossing() function: seems like this is making first BC generation, but neither Material has cross symbol')
        raise ValueError('Inside backCrossing() function: seems like this is making first BC generation, but neither Material has cross symbol')
     
            
    #Case of later generation backcross of inbred female by BCx male
    #seed_gen passed from seed_generation() is not correct in this case, need to update it
    elif bc_gen == 0 and match_bc: 
        print('made it to later gen bc of inbred female by BCx male')
        update_bc_gen = int(pair_bc_gen) + 1
        recurrent_doses = '*' + str(update_bc_gen + 1)
        seed_gen = "BC" + str(update_bc_gen) + "F1"
        Material = re.sub('\*[0-9]+', recurrent_doses, pair['Material'][0]) #eg replace (B73*2)Oh43 with (B73*3)Oh43

    #Case of later generation backcross of BCx female by inbred male
    elif bc_gen != 0:
        print('backCrossing case 4')
        update_bc_gen = int(bc_gen) + 1
        recurrent_doses = '*' + str(update_bc_gen + 1)
        Material = re.sub('\*[0-9]+', recurrent_doses, r['Material']) #eg replace (B73*2)Oh43 with (B73*3)Oh43
        
    #print('Material:')
    #print(Material)    
    
    seed_ID = r['Row_ID'] + u'\xd7' + pair['Row_ID'][0][4:]
    source_ID = '(' + r['Source ID'] + ')' + u'\xd7' + '(' + pair['Source ID'] + ')'
      
    #Bulk shelling
    if r['shell'] == 'BULK': 
        seed_ID = seed_ID + '-BLK'
        return DataFrame({'seed_ID' : seed_ID,
                          'Material' : Material,
                          'seed_gen' : seed_gen,
                          'ear_number': earno_cross,
                          'ear_quality' : earq_cross,
                          'source_ID' : source_ID}, index = [0])

    else:
        result = DataFrame()
        for ear in range(1, int(earno_cross) + 1):           
            ear_ID = seed_ID + '-' + str(ear)
            ear_to_row = DataFrame({'seed_ID' : ear_ID,
                          'Material' : Material,
                          'seed_gen' : seed_gen,
                          'ear_quality' : earq_cross,
                          'source_ID' : source_ID}, index = [0])
            result = result.append(ear_to_row)

        return result
    
    print("at end of backCrossing() with no return value")
        #
        
#Function to generate labels for free crossing. 
#updated so that free_pairs is a nested dict with this structure: {female:male} with male as a nested dict: {male:(earno, earq)}
#We require that female and male keys are type int, this will fail otherwise
def freeCrossing(r, inputdf, seed_gen):
    row = r['Row'] #this will be int type
    result = DataFrame()
    print('inside freeCrossing(), row = ')
    print(row)
                                   
    for pair_row in free_pairs[row]: #each pair_row is the key for a nested dict of {male:(earno, earq)}
        pair = inputdf[inputdf['Row'] == pair_row] # pull the paired row as a Series from the inputdf data frame
        pair.index = [0] #set the index to zero so we can easily access the contents  
        print('inside freeCrossing(), pair = ')
        print(pair)
        if r['Material'] == pair['Material'][0]: 
            Material = r['Material']
        else: 
            Material = parenth(r['Material']) + u'\xd7' + parenth(pair['Material'][0]) #access the 0 index to get just the Material content itself instead of the Series
        seed_ID = r['Row_ID'] + u'\xd7' + pair['Row_ID'][0][4:]
        source_ID = '(' + r['Source ID'] + ')' + u'\xd7' + '(' + pair['Source ID'] + ')'
        earno_cross = free_pairs[row][pair_row][0]
        earq_cross = free_pairs[row][pair_row][1]
        
        #Bulk shelling
        if r['shell'] == 'BULK': 
            seed_ID = seed_ID + '-BLK'
            bulk_row = DataFrame({'seed_ID' : seed_ID,
                              'Material' : Material,
                              'seed_gen' : seed_gen,
                              'ear_number': earno_cross,
                              'ear_quality' : earq_cross,
                              'source_ID' : source_ID}, index = [0])
            result = result.append(bulk_row)

        else:

            for ear in range(1, int(earno_cross) + 1):
                ear_ID = seed_ID + '-' + str(ear)
                ear_to_row = DataFrame({'seed_ID' : ear_ID,
                              'Material' : Material,
                              'seed_gen' : seed_gen,
                              'ear_quality' : earq_cross,
                              'source_ID' : source_ID}, index = [0])
                result = result.append(ear_to_row)
    return result

def intermating(r, seed_gen, earno_cross, earq_cross):
    print('inside intermating for row ' + str(r['Row']))

    #Bulk shelling
    if r['shell'] == 'BULK':  
        #Generate seed_ID for case of bulk shelling
        seed_ID = r['Row_ID'] + '-BLK'
        return DataFrame({'seed_ID' : seed_ID,
                          'Material' : r['Material'],
                          'seed_gen' : seed_gen,
                          'ear_number': earno_cross,
                          'ear_quality' : earq_cross,
                          'source_ID' : r['Source ID']}, index = [0])
 
    
    #Single ear shelling        
    else: 
        result = DataFrame()
        for ear in range(1, int(earno_cross) + 1):
            seed_ID = r['Row_ID'] + '-' + str(ear)
            ear_to_row = DataFrame({'seed_ID' : seed_ID,
                                    'Material' : r['Material'],
                                    'seed_gen' : seed_gen,
                                    'ear_quality' : earq_cross,
                                    'source_ID' : r['Source ID']}, index = [0])
            result = result.append(ear_to_row)
        return result

    
#Function to generate self and cross generations from same rows.
def selfAndCross(r, inputdf, bc_gen, rm_gen, prefix_char, prefix_num, suffix, earno_self, earq_self, earno_cross, earq_cross, cross_type, seed_gen = 0):
    #the trick here will be to make two output rows for each input row, one for selfing, other for crossing
    #bind them together into two rows of a data frame that is returned
    selfs = selfing(r, bc_gen, rm_gen, prefix_char, prefix_num, suffix, earno_self, earq_self)
    if cross_type == "CR": cross = crossing(r, inputdf, 'F1', earno_cross, earq_cross)
    elif cross_type == "BC": cross = backCrossing(r, inputdf, bc_gen, seed_gen, earno_cross, earq_cross)
    elif cross_type == "FC": cross = freeCrossing(r, inputdf, "F1")
    both = pd.concat([selfs, cross], axis = 0, ignore_index=True, sort = True)
    return both
    
#Function to generate sib and cross generation from same rows.
def sibAndCross(r, inputdf, bc_gen, rm_gen, prefix_char, prefix_num, suffix, sib, earno_self, earq_self, earno_cross, earq_cross, cross_type):
    #the trick here will be to make two output rows for each input row, one for selfing, other for crossing
    #bind them together into two rows of a data frame that is returned
    sibbed =  sibbing(r, bc_gen, rm_gen, prefix_char, prefix_num, suffix, sib, earno_self, earq_self)
    if cross_type == "CR": cross = crossing(r, inputdf, 'F1', earno_cross, earq_cross)
    elif cross_type == "BC": cross = backCrossing(r, inputdf, bc_gen, seed_gen, earno_cross, earq_cross)
    elif cross_type == "FC": cross = freeCrossing(r, inputdf, seed_gen)
    both = pd.concat([sibbed, cross], axis = 0, ignore_index=True, sort = True)
    return both
    

#Function to take decomposed plot generation and row index and make the seed generation:
def seed_generation(r, inputdf, bc_gen, rm_gen, prefix_char, prefix_num, suffix, sib):
    '''
    Takes the parsed plot generation information from split_generation function and passes to a function specific for the pollination type
    '''
    print('inside seed_generation')
    print(r)

    #print('inside seed_generation for row ' + str(r['Row']))
    #print(r)

    earno_self = r['earno_self'] 
    earq_self = r['earq_self'] 
    earno_cross = r['earno_cross'] 
    earq_cross = r['earq_cross'] 
    if r['Poll_Type'] == "SF" and r['earno_self'] > 0: 
        return selfing(r, bc_gen, rm_gen, prefix_char, prefix_num, suffix, earno_self, earq_self) #no need to pass sibbing generation if selfing, F6:7Sib2 selfed is same as F6:7 selfed
    elif r['Poll_Type'] == "SB" and r['earno_self'] > 0: 
        return sibbing(r, bc_gen, rm_gen, prefix_char, prefix_num, suffix, sib, earno_self, earq_self)

    #Crossing always returns "F1" generation. We now separately specify TC for topcrossing and BC for backcrossing, 
    #IN for intermating blocks of rows, and FC for free crosses between non-adjacent rows
    elif r['Poll_Type'] == "CR" and r['earno_cross'] > 0:
        seed_gen = "F1"
        return crossing(r, inputdf, seed_gen, earno_cross, earq_cross)
    elif r['Poll_Type'] == "TC" and r['earno_cross'] > 0:
        seed_gen = "TC"
        return crossing(r, inputdf, seed_gen, earno_cross, earq_cross)
    elif r['Poll_Type'] == "BC" and r['earno_cross'] > 0: 
        seed_gen = 'BC' + str(bc_gen + 1) + 'F1'
        return backCrossing(r, inputdf, bc_gen, seed_gen, earno_cross, earq_cross)
    elif r['Poll_Type'] == "FC" and r['earno_cross'] > 0: 
        seed_gen  = "F1"
        return freeCrossing(r, inputdf, seed_gen)
    elif r['Poll_Type'] == "IN" and r['earno_cross'] > 0: 
        seed_gen = 'RM' + str(rm_gen + 1) + 'S0'
        return intermating(r, seed_gen, earno_cross, earq_cross) 
    elif r['Poll_Type'] == "SF/CR" or r['Poll_Type'] == "CR/SF" : 
        if r['earno_self'] > 0 and r['earno_cross'] == 0: 
            return selfing(r, bc_gen, rm_gen, prefix_char, prefix_num, suffix, earno_self, earq_self) 
        elif r['earno_self'] == 0 and r['earno_cross'] > 0:
            seed_gen = "F1"
            return crossing(r, inputdf, seed_gen, earno_cross, earq_cross)
        else:
            return selfAndCross(r, inputdf, bc_gen, rm_gen, prefix_char, prefix_num, suffix,earno_self, earq_self, earno_cross, earq_cross, "CR") 
    elif r['Poll_Type'] == "SB/CR" or r['Poll_Type'] == "CR/SB" : 
        if r['earno_self'] > 0 and r['earno_cross'] == 0: 
            return  sibbing(r, bc_gen, rm_gen, prefix_char, prefix_num, suffix, sib, earno_self, earq_self) 
        elif r['earno_self'] == 0 and r['earno_cross'] > 0:
            seed_gen = "F1"
            return crossing(r, inputdf, seed_gen, earno_cross, earq_cross)
        else:
            return sibAndCross(r, inputdf, bc_gen, rm_gen, prefix_char, prefix_num, suffix, sib, earno_self, earq_self, earno_cross, earq_cross, "CR") 
    elif r['Poll_Type'] == "SF/BC" or r['Poll_Type'] == "BC/SF" : 
        if r['earno_self'] > 0 and r['earno_cross'] == 0: 
            return selfing(r, bc_gen, rm_gen, prefix_char, prefix_num, suffix, earno_self, earq_self) 
        elif r['earno_cross'] > 0:
            seed_gen = 'BC' + str(bc_gen + 1) + 'F1'
            if r['earno_self'] == 0:            
                return backCrossing(r, inputdf, bc_gen, seed_gen, earno_cross, earq_cross)
            else:
                return selfAndCross(r, inputdf, bc_gen, rm_gen, prefix_char, prefix_num, suffix,earno_self, earq_self, earno_cross, earq_cross, "BC", seed_gen) 
    elif r['Poll_Type'] == "SF/FC" or r['Poll_Type'] == "FC/SF" : 
        if r['earno_self'] > 0 and r['earno_cross'] == 0: 
            return selfing(r, bc_gen, rm_gen, prefix_char, prefix_num, suffix, earno_self, earq_self) 
        elif r['earno_self'] == 0 and r['earno_cross'] > 0:
            seed_gen = "F1"
            return freeCrossing(r, inputdf, seed_gen, earno_cross, earq_cross)
        else:
            return selfAndCross(r, inputdf, bc_gen, rm_gen, prefix_char, prefix_num, suffix,earno_self, earq_self, earno_cross, earq_cross, "FC") 

    print("at end of seed_generation with no return value")
    
#Make a function to decompose generation into 8 parts potentially:
#BCwFx:ySbz
#BC = BC: designator for a backcross generation
#w = BC_gen: number of backcross generations
#F = prefix_char: prefix character (F, S, or RM for random-mated)
#x = prefix_num: generation number of last common ancestor
#y = suffix: current generation
#Sb = Sib: designator for a sibbing generation for seed maintenance (not inbreeding)
#z = sib: number of sib generations

def split_generation(r, inputdf):
    plot_gen = r['Gen'] # this gets the current generation of a plot
    print('inside split_generation for row ' + str(r['Row']))
    print('plot_gen is ' + str(plot_gen))
    match_bc = re.match('BC[0-9]', plot_gen) #returns match object only if plot_gen starts with 'BC'
    #if there was a BC generation designator, we now get the BC generation number, then strip off the BC  prefix and deal with the rest
    if match_bc:
        bc_gen = int(plot_gen[2:match_bc.end()]) #this gets the bc generation number 
        plot_gen = plot_gen[match_bc.end():] #this is everything after BC
    else: bc_gen = 0
 
    match_rm = re.match('RM[0-9]', plot_gen) #returns match object only if plot_gen starts with 'RM'
    #if there was a RM generation designator, we now get the RM generation number, then strip off the RM  prefix and deal with the rest
    if match_rm:
        rm_gen = int(plot_gen[2:match_rm.end()]) #this gets the bc generation number 
        plot_gen = plot_gen[match_rm.end():] #this is everything after RM
    else: rm_gen = 0
   
    split_gen = re.split(':', plot_gen) #this will split plot_gen into a list of two if ':' in generation, else is a list of one
    match_pre = re.match('[a-zA-Z]+',split_gen[0]) #match the prefix letters
    
    if match_pre: 
        prefix_char = match_pre.group() #this gets the prefix characters from the match object
        prefix_num = (split_gen[0][match_pre.end():]) #this pulls everything after the characters before the :, which should be numbers of prefix gen    

    #catch cases where the prefix character is not alphabetical, like ? or something weird
    else: 
        prefix_char = '?'
        prefix_num = 0

    #this is a little tricky, sometimes you have a prefix num (need to convert to int), other times not (converting to int raises error)
    if prefix_num: prefix_num = int(prefix_num)
    else: prefix_num = 0

    if len(split_gen) > 1:
        suffix = split_gen[1]
        #check if 'sib' is appended to suffix:
        match_sib = re.search('Si?b',suffix) # will match 'Sib' or 'Sb'
        
        #now partition suffix into the suffix number and sib number, and also convert the numbers to ints if they exist
        if match_sib: 
            suffix_num = suffix[:match_sib.start()]
            sib = int(suffix[match_sib.end():])
        else: 
            suffix_num = suffix
            sib = 0
        if len(suffix_num) > 0 : suffix_num = int(suffix_num)
        return seed_generation(r, inputdf, bc_gen, rm_gen, prefix_char, prefix_num, suffix_num, sib)
    else:
        return seed_generation(r, inputdf, bc_gen, rm_gen, prefix_char, prefix_num, suffix = 0, sib = 0)
        
#We need a function labeller that will accept as an argument a single row of the input data frame and return a new data frame with potentially multiple rows (because the several ears harvested from one row will be come several rows in the harvest shelling label data frame).
#The labeller function serves as the top level of determining how each row of the input data frame will be processed. In most cases, it will just pass the information to 'split_generation', but it serves to do some checking of the input information and handling unusual cases.
#Unfortunately, we cannot simply use DataFrame.apply() method to apply this function to each row of the data frame. The reason is that DataFrame.apply() will apply a function to each row, but it expects as a return type a Series representing a single row of same dimensions as the input row. This makes sense when you are using .apply() to make a transformation of the input data frame row by row. But instead we are making a new object that will have a different number of rows (and a different number of columns) than the input data frame.
#The simplest way to do this is to create a groupby object from the original data frame. We can group by individual rows. Then we can use .apply() method on the groupby object, because this method allows us to return different dimensions than the input object. It will also concatenate all of the groupby() return objects into a single data frame, this is what we want.
#Note a tricky point here, if we pass a groupby object directly to split_generation, we are passing a data frame of the groupby object with one row. Split_generation() is designed to accept a Series not a Data Frame. The first line inside of split_generation is: plot_gen = r['Gen'] . If we pass a DataFrame as argument r to split-generation, this line produces a SERIES not a string, so something like: 0 'F6:7' instead of 'F6:7'. And that causes the function to fail.
#To be sure we get it right, we first take the input groupby object 'group' and turn it into a Series called 'r'. This is done at beginning of labeller function with code r = group.iloc[0]. I think r = group.irow(0) would give same result. Then we pass the Series r as the input argument to split_generation(r).
#In addition, it seems to be that if the first row of the original data frame has a zero count for ears, then the labels df is not built up properly. So, we need to first scan the group and if the ear counts are both zero, return a data frame with placeholders.

def labeller(group, inputdf):
    r = group.iloc[0]
    #r is a Series that represents a single row of the data frame
    print('labeller for row ' + str(r['Row'])) #useful for debugging
    #check if both ear counts are zero
    
    if r['earno_self'] == 0 and r['earno_cross'] == 0:
        return DataFrame({'seed_ID' : 'NO',
                          'Material' : 'NO',
                          'seed_gen' : 'NO',
                          'ear_quality' : 'NO',
                          'source_ID' : 'NO', 
                          'shell' : 'NO'}, index = [0])

    elif r['shell'] == 'MULTI-BULK':
        #switch shell to 'BULK' to be able to pass to split_generation function
        r['shell'] = 'BULK'
        label = split_generation(r, inputdf)
        if not label.empty: 
            label['shell'] = r['shell']
            return label        
        else: #use this to catch any problem rows
            print('split_generation returned None at row:')
            print(r['Row'])
            print(r)
    
    else:
        label = split_generation(r, inputdf)
        print("current label = " + str(label))
        if not label.empty: 
            label['shell'] = r['shell']
            return label        
        else: #use this to catch any problem rows
            print('split_generation returned None at row:')
            print(r['Row'])
            print(r)
   

#createLabels function will take the input file, do some standard processing, group it by rows, then pass the groupby object to the labeller function. It then will do a little processing on the output to split into appropriate sub-dataFrames and write those to hard drive.
def createLabels(inputdf, paired_rows, nursery_prefix, free_pairs_dict = {}, multirows = {}, multis_list = {}, return_labels = False):

    global nursery_pre
    nursery_pre = nursery_prefix
    global free_pairs
    free_pairs = free_pairs_dict
    
    #Create a global dictionary called 'pairs' that will be accessed by some of the label generating functions to identify the paired row companions for crossing and backcrossing.
    global pairs
    pairs = dict()
    pair = 0
    for i in paired_rows:
        start, stop = i
        for row in range(start, stop + 1, 2):
            pairs[row] = pair%2
            pairs[row + 1] = pair%2
            pair += 1
    
    inputdf['Row_ID'] = inputdf.apply(rowToRowID, axis = 1)
    
    #Need to convert Row to integer type
    inputdf[['Row']] = inputdf[['Row']].astype(int)
    
    #Need to correct/change some Poll_Type values. For consistency, 'CR/SF' should always be 'SF/CR'. 
    #Implement global replacements with replace method of data frames.
    inputdf['Poll_Type'].replace('CR/SF', 'SF/CR', inplace = True)
                
#Create a global dictionary called 'pairs' that will be accessed by some of the label generating functions to identify 
#the paired row companions for crossing and backcrossing.
    pairs = dict()
    pair = 0
    for i in paired_rows:
        start, stop = i
        for row in range(start, stop + 1, 2):
            pairs[row] = pair%2
            pairs[row + 1] = pair%2
            pair += 1
    
#Partition the data frame into two pieces, one for individual or paired row shelling, and a second data frame for multi-row bulks. They will require different approaches to generating the labels.    
    singles_pairs = inputdf.loc[(inputdf['earno_self'] != 0) | (inputdf['earno_cross'] != 0)]
    
#Convert the data frame to a groupby object, grouping by rows. This will allow us to use the apply method to the groupby object 
#to make a new data frame containing the shelling labels.
    counts_by_row = singles_pairs.groupby('Row', group_keys=False, as_index=False)
    
#Apply the labeller function to each groupby object (row of counts_by_rows data frame)
    labels = counts_by_row.apply(labeller, inputdf) 
    
#post-processing
    labels['Row'] = labels['seed_ID'].map(lambda x: int(x[4:8]))
    
#Partition the data frame into single ear and bulk shelling lists
    singles = labels.loc[labels['shell'] == 'SINGLE-EAR']
    bulks = labels.loc[labels['shell'] == 'BULK']
    bulks = bulks.loc[-bulks['Row'].isin(multis_list)]

#Now take out the multi-bulk rows and combine them into consecutive multi-row bulks according to the multi-rows list
    multis = labels.loc[labels['Row'].isin(multis_list)]

    if not multis.empty:
        for tup in multirows:
            print('Current start/end pair in multirows:')
            print(tup)
            rowset = multis[multis['Row'].isin(range(tup[0], tup[1]+1))]
            rowset0 = rowset.iloc[[0]] #take only the first row
            rowset0['seed_ID'] = rowset0.loc[rowset.index[0],'seed_ID'][0:8] + '-' + str(tup[1]).zfill(4) + '-BLK'
            #check if multis2 exists, if it does concat, if not (first element in multirows), then create it
            try:
                multis2 = pd.concat([multis2, rowset0], axis = 0, ignore_index = True, sort = True)
            except NameError:
                multis2 = rowset0
                
#Write the labels data frame to a csv file. Automatically attach current date to file name (provides some protection against 
#writing over current files on hard drive).
    date = strftime("%m-%d-%Y") #this gets current date as MM-DD-YYYY string

    try:
        if not singles.empty: 
            filenameS = nursery_pre + ' Single Ear Shelling Labels ' + date + '.csv' #filename for single ear shelling list
            singles.to_csv(filenameS, index = False, index_label=False, columns = ['seed_ID', 'Material', 'seed_gen', 'source_ID', 'ear_quality', 'shell'], encoding = 'cp1252')

        if not bulks.empty:    
            filenameB = nursery_pre + ' Bulk Shelling Labels ' + date + '.csv' #filename for bulk shelling list
            bulks.to_csv(filenameB, index = False, index_label=False, columns = ['seed_ID', 'Material', 'seed_gen', 'source_ID', 'ear_number', 'ear_quality', 'shell'], encoding = 'cp1252')
        
        if not multis.empty:
            filenameM = nursery_pre + ' Multirow Shelling Labels ' + date + '.csv' #filename for multirow shelling list
            multis2.to_csv(filenameM, index = False, index_label=False, columns = ['seed_ID', 'Material', 'seed_gen', 'source_ID', 'ear_number', 'ear_quality', 'shell'])

    except NameError:
        pass #just silently pass over exception that occurs if one of these data frames does not exist
    
    if return_labels:
        return labels

#########################################################################################################
