from molmass import Formula, FormulaError
import os
import re
import fitz
from fpdf import FPDF

import sys
import math
import json
import string
import pandas as pd
from datetime import date

def calculateError(found_mass_from_si, calculated_mass):
    return abs(round((calculated_mass / found_mass_from_si - 1) * 10 ** 6, 1))

def setErrorLevel(error_level, new_level):
    if error_level < new_level:
        return new_level   
    else:
        return error_level

# SI -- calculated and reported masses differ in their integers
# Input:
#     - calculated mass from si
#     - found mass from si
# Return: 
#     - True = there is a different integer
#     - False = otherwise    
def checkDifferentIntegers(calculated_mass_from_si, found_mass_from_si):

    calculated_mass_from_si_str = str(int(calculated_mass_from_si))
    found_mass_from_si_str = str(int(found_mass_from_si))

    typo_si = False
    for char1, char2 in zip(calculated_mass_from_si_str, found_mass_from_si_str):
        if char1 == '.' or char2 == '.':
            break
        if char1 != char2:
            return True
    return False

# Input:
#     - calculated mass from si
#     - found mass from si
# Return: 
#     - True = there are swapped integers excluding the last two digits
#     - False = otherwise    

def checkSwappedIntegers(calculated_mass_from_si, found_mass_from_si):
    # Convert both numbers to strings (excluding the last character)
    str_num1 = str(calculated_mass_from_si)[:-1]
    str_num2 = str(found_mass_from_si)[:-1]

    # If lengths of the numbers are different, return False
    if len(str_num1) != len(str_num2):
        return False

    # Find mismatched positions and check for a valid swap
    mismatches = []
    for i in range(len(str_num1)):
        if str_num1[i] != str_num2[i]:
            mismatches.append(i)
            # If there are exactly 2 mismatches, check for a valid swap
            if len(mismatches) == 2:
                i, j = mismatches
                if str_num1[i] == str_num2[j] and str_num1[j] == str_num2[i]:
                    return True
            # If more than 2 mismatches, it's not a valid swap
            elif len(mismatches) > 2:
                return False

    # If we finish the loop without finding exactly 2 mismatches, return False
    return False


# typo in calculated and/or measurement / either differ from the calculated mass by one digit
# last character must be the same or will reject
# Input:
#     - calculated mass from si
#     - found mass from si
#     - calculated mass from neutral/cation/anion
# Return: 
#     - 0 = not a typo problem
#     - 1 = typo with calculated from si only
#     - 2 = typo with found from si only
#     - 3 = typo with both calculated from si and found from si   
def checkTypo(calculated_mass_from_si, found_mass_from_si, calculated_mass):
    calculated_mass_from_si_str = "{:.4f}".format(calculated_mass_from_si)[-4:-2]
    found_mass_from_si_str = "{:.4f}".format(found_mass_from_si)[-4:-2]
    calculated_mass_str = "{:.4f}".format(calculated_mass)[-4:-2]

    typo_si = False
    for char1, char2 in zip(calculated_mass_from_si_str, found_mass_from_si_str):
        if char1 != char2:
            if typo_si:
                typo_si = False
                break
            else:
                typo_si = True
    
    typo_calculated = False

    return typo_si or typo_calculated

# check if the mass of the added ion was wrongly added to the neutral molecule
# output:
#   0 = not because of the ion + neutral mass
#   1 = neutral + assumed mass = calculated mass from si
#   2 = neutral - actual mass + assumed mass = calculated mass from si
def checkIon(molecular_formula, calculated_mass_from_si, calculated_mass_from_neutral, added_ion, assumed_mass):
    # check if the ion is in the molecular formula
    elements = Formula(added_ion).composition()
    composition = Formula(molecular_formula).composition()
    for curr_element, content in elements.items():
        if curr_element != 'e-':
            if curr_element not in composition:
                return 0
            if content.count < elements[curr_element].count:
                return 0
    # neutral + assumed mass
    if calculated_mass_from_neutral + assumed_mass == calculated_mass_from_si:
        return 1
    
    # neutral - actual mass of added ion + assumed mass    
    if float("{:.4f}".format(calculated_mass_from_neutral - Formula(added_ion).monoisotopic_mass + assumed_mass)) == calculated_mass_from_si:
        return 2

# Input:
#     - molecular formula
#     - calculated mass from si
#     - found mass from si
#     - element to add
#     - upper range for how many of 'element to add' to add
# Return: 
#     - number of element(s) to add s.t. the new mass with the addition == calculated mass from si
def checkAddition(molecular_formula, calculated_mass_from_si, element_to_add, upper_range):
    calculated_mass = Formula(molecular_formula).monoisotopic_mass

    for i in range(upper_range+1):
        new_calculated_mass = float("{:.4f}".format(calculated_mass + (i) * Formula(element_to_add).monoisotopic_mass))
        if new_calculated_mass == calculated_mass_from_si:
            return i

    return 0

# Input:
#     - molecular formula
#     - calculated mass from si
#     - found mass from si
#     - element to remove
#     - upper range for how many of 'element to remove' to remove
# Return: 
#     - number of element(s) to remove s.t. the new mass with the addition == calculated mass from si
def checkRemove(molecular_formula, calculated_mass_from_si, element_to_remove, upper_range):
    calculated_mass = Formula(molecular_formula).monoisotopic_mass
    
    # check if element_to_remove is in molecular formula. element_to_remove could be 1 or more elements
    elements = Formula(element_to_remove).composition()
    composition = Formula(molecular_formula).composition()
    for curr_element, content in elements.items():
        if curr_element != 'e-':
            if curr_element not in composition:
                return 0
            else:
                if composition[curr_element].count < content.count * upper_range:
                    upper_range = math.floor(composition[curr_element].count / content.count)

    for i in range(upper_range+1):
        new_calculated_mass = float("{:.4f}".format(calculated_mass - (i) * Formula(element_to_remove).monoisotopic_mass))
        if new_calculated_mass == calculated_mass_from_si:
            return i

    return 0

def checkReplace(molecular_formula, calculated_mass_from_si, element_to_remove, element_to_add):
    calculated_mass = Formula(molecular_formula).monoisotopic_mass
    composition = Formula(molecular_formula).composition()

    if element_to_remove in composition:
        new_calculated_mass = float("{:.4f}".format(calculated_mass - Formula(element_to_remove).monoisotopic_mass + Formula(element_to_add).monoisotopic_mass))
        return new_calculated_mass == calculated_mass_from_si
    else:
        return False

# composition -> molecular formula
def compositionToFormula(composition, measuring_mode):
    formula_str = []
    for element in composition:
        if element != 'e-':
            if composition[element]['count'] != 0:
                if composition[element]['count'] == 1:
                    formula_str.append(f'{element}')
                else:
                    formula_str.append(f"{element}{composition[element]['count']}")
    molecular_formula_neutral = '['+''.join(formula_str)+']'
    if measuring_mode == "cation":
        return molecular_formula_neutral+"+"
    elif measuring_mode == "anion":
        return molecular_formula_neutral+"-"
    else:
        return molecular_formula_neutral

def compositionToDict(composition):
    return {element: {"count": content.count} for element, content in composition.items()}


# when you run the program, python 'filepath' threshold number_of_files output_filepath fiilepaths ...
threshold = float(sys.argv[1])
numfiles = int(sys.argv[2])
outfile = sys.argv[3]
filepaths = [None] * numfiles
neutral = sys.argv[4]
if neutral == '0':
    neutral = False
elif neutral == '1':
    neutral = True
for i in range(numfiles):
    filepaths[i] = sys.argv[i+5]

file_index = 0
outputJSONFile = open(outfile, 'w')

for filepath in filepaths:
    title = "title"
    file_index = file_index + 1
    file_request = []

    total_examples = 0
    incorrect_examples = [0, 0, 0, 0] # indices: 0 = worst, 1 = fixable, 2 = minor
    invalid_inputs = 0
        
    try:
        # extract text from pdf
        if os.path.basename(filepath).lower() == 'desktop.ini':
            file_contents_arr = []
        else:
            pdf_document = fitz.open(filepath)
            text_content = []

            for page_num in range(pdf_document.page_count):
                page = pdf_document.load_page(page_num)
                text_content.append(page.get_text()) # index = 0 => page = 1
            file_contents_arr = text_content
        
        # pandas df to hold the extracted data
        extracted_data = pd.DataFrame(columns=['page number', 'extracted text', 'initial comment', 'molecular ion type',
                                               'molecular formula neutral', 'molecular formula cation', 'molecular formula anion', 
                                               'measuring mode', 'sodium', 'calculated mass from si', 'found mass from si', 
                                               'calculated mass from neutral', 'calculated mass from cation', 'calculated mass from anion', 
                                               'mass error from si', 'mass error from neutral', 'mass error from cation', 'mass error from anion'])
        
        # line structure to determine where in the line the necessary words are
        line_structure = { # default structure -- change as needed
            "molecular formula": 2,
            "calculated mass from si": 3,
            "found mass from si": 5,
            "fixed": False
        }

        curr_string = "" # to help with any carry-over from the previous page
        
        for page_num in range(0, len(file_contents_arr)):
            print("page", page_num)
            file_contents = curr_string.rstrip() + file_contents_arr[page_num].lstrip() 
            # print(file_contents)
            try:
                # check if this page contains any hrms data
                file_contents_removed = re.sub(r'[^a-zA-Z0-9]', '', file_contents)

                if not 'hrms' in file_contents_removed.lower() and not 'calc' in file_contents_removed.lower() and not 'found' in file_contents_removed.lower():
                    raise ValueError('hrms not found on page') # will throw valueerror if not and continue

                # remove short lines
                lines = file_contents.splitlines()
                lines = [line if len(line) >= 7 else '\n' for line in lines]
                file_contents = '\n'.join(lines)

                if curr_string == "HRMS":
                    print("HRMS curstring")
                    print(file_contents)
                
                # search for all hrms data in the page
                curr_index = -1
                hrms_index = -1
                while curr_index < len(file_contents):
                    # print(curr_index)
                    # print(curr_index)
                    # print(file_contents[curr_index+1:])
                    hrms_index = file_contents.lower().find("hrms", curr_index+1)
                    curr_index = file_contents.lower().find("cal", curr_index+1)
                    # print(hrms_index)
                    print(curr_index, hrms_index)
                    print(file_contents[curr_index:curr_index+100])
                    
                    if curr_index == -1:
                        if hrms_index != -1:
                            # print("hrms index", hrms_index)
                            new_line_index = file_contents.find('\n', hrms_index)
                            # print()
                            curr_string = file_contents[hrms_index:new_line_index]
                        
                        print("breaking ?")
                        break # no more hrms data left on the page
                    # elif curr_index == -1 and hrms_index != -1:
                        # curr_string = file_contents[len(file_contents)-100:len(file_contents)]
                        # print(curr_string)
                        # curr_string = "hrms "
                        # break
                        # curr_index = curr_index + 100
                    else:
                        try:
                            # check that this line is actually hrms data
                            start_index = curr_index - 50
                            end_index = curr_index + 100
                            if curr_index-50 < 0:
                                start_index = 0
                            if curr_index+100 >= len(file_contents):
                                end_index = len(file_contents)
                            found_string_test = file_contents[start_index:end_index]
                            found_string_removed = re.sub(r'[^a-zA-Z0-9]*', '', found_string_test).lower() 

                            # print(found_string_test)
                            # print(found_string_removed) 
                            
                            if ('hrms' in found_string_removed and 'found' in found_string_removed):
                                print("hrms in found string")
                                if curr_string != "":
                                    curr_string = ""

                                found_string_removed = re.sub(r'[^a-zA-Z0-9+-\[\]\(\)]*', '', file_contents[start_index:end_index]).lower() 
                                print(found_string_removed)

                                # determine whether cation/anion/neutral measuring mode
                                measuring_mode = ""
                                molecular_ion_type = 'unknown'
                                sodium = False

                                cation_mit_index = found_string_removed.find(']+', 0)
                                if cation_mit_index == -1:
                                    found_string_removed.find(')+', 0)

                                anion_mit_index = found_string_removed.find(']-', 0)
                                if anion_mit_index == -1:
                                    found_string_removed.find(')-', 0)

                                if cation_mit_index != -1 and ((cation_mit_index < anion_mit_index and anion_mit_index != -1) or (anion_mit_index == -1)): 
                                    molecular_ion_type = 'cation'
                                elif cation_mit_index != -1 and ((anion_mit_index < cation_mit_index and cation_mit_index != -1) or (cation_mit_index == -1)):
                                    molecular_ion_type = 'anion'
                                else:
                                    molecular_ion_type = 'unknown'

                                # m+h -> cation
                                # m+ -> cation
                                # m- -> anion
                                # m-h -> anion

                                # cation
                                # m+nh4
                                # m+Na
                                # m+k

                                # anion
                                # m+cl
                                # m+ch3coo

                                # clean up the found_string
                                # found_string = re.sub('\n+', '\n', file_contents[curr_index:curr_index+75])
                                found_string = re.sub(r'\s+', ' ', found_string_test)
                                # print(found_string)
                                found_string = re.sub(r'[:;,]', ' ', found_string)
                                found_string = re.sub(r'[\+-]+', '', found_string)
                                found_string = re.sub(r'\(.*?\)', '', found_string)
                                found_string = re.sub(r'\[.*?\]', '', found_string)
                                # print(found_string)
                                found_string = re.sub(r'h\s*r\s*m\s*s', 'hrms', found_string, flags=re.IGNORECASE)

                                found_string = found_string.replace('[M+Na]+', '')
                                found_string = found_string.replace('Na]+', '')
                                found_string = found_string.replace('(M+Na)+', '')
                                found_string = found_string.replace('(M + Na)+', '')
                                found_string = found_string.replace('[M+H]+', '')
                                found_string = found_string.replace('[M-H]+', '')
                                found_string = found_string.replace('[M+H]', '')
                                found_string = found_string.replace('[M]', '')
                                found_string = found_string.replace('[M + H]+', '')
                                found_string = found_string.replace('[M - H]+', '')
                                found_string = found_string.replace('m/z', '')
                                found_string = found_string.replace('m/s', '')
                                found_string = found_string.replace('(M+H)+', '')
                                found_string = found_string.replace('(M + H)+', '')
                                found_string = found_string.replace('[]', '')
                                found_string = found_string.replace('[M', '')
                                found_string = found_string.replace('H]+', '')

                                found_string = re.sub(r'(\s+)(\d+).(\s+)(\d+)(\s+)', r'\1\2.\4\5', found_string)
                                found_string = re.sub(r'\d+[\+\-–]', '', found_string)

                                # tokenize the found string
                                str_split = re.split(r'\s+', found_string)
                                str_split = [curr_str for curr_str in str_split if len(curr_str) > 1]

                                # print(str_split)

                                for i in range(len(str_split)):
                                    if 'hrms' in re.sub(r'[^a-zA-Z0-9]*', '', str_split[i]).lower():
                                        # print("hrms found at", i)
                                        str_split = str_split[i:]
                                        str_split[0] = 'hrms'
                                        break
                                
                                # print(str_split)

                                # establish line_structure for the first instance to find where 
                                # molecular formula, calculated mass from si, and found mass from si are in the tokenized string
                                between_hrms_found = 0
                                if not line_structure['fixed']:
                                    for i in range(len(str_split)):
                                        # print(i, str_split[i])
                                        if between_hrms_found == 0 and 'hrms' in re.sub(r'[^a-zA-Z0-9]*', '', str_split[i]).lower():
                                            between_hrms_found = between_hrms_found + 1
                                        elif between_hrms_found == 1 and str_split[i].lower() == 'found':
                                            between_hrms_found = between_hrms_found + 1
                                            line_structure.update({'found mass from si': i+1})
                                        if between_hrms_found == 1:
                                            try:
                                                calculated_mass_from_si = str_split[i]
                                                calculated_mass_from_si = calculated_mass_from_si.rstrip('.')
                                                # print(calculated_mass_from_si)
                                                calculated_mass_from_si = float(calculated_mass_from_si)

                                                line_structure.update({'calculated mass from si': i})
                                            except:
                                                # continue
                                                pass
                                            
                                            try:
                                                mass = float("{:.4f}".format(Formula(str_split[i]).monoisotopic_mass))
                                                line_structure.update({'molecular formula': i})
                                            except:
                                                # continue     
                                                pass
                                    line_structure.update({'fixed': True})
                                    print(line_structure)
                                if len(str_split) < line_structure['found mass from si']: # if the found_string is cut off, continue to next page
                                    curr_string = found_string
                                elif not re.search(r'\d', str_split[0]):
                                    print(str_split)

                                    molecular_formula = re.sub(r'\W+', '', str_split[line_structure['molecular formula']]).rstrip('+-.[]')
                                    molecular_formula_cation = '[' + molecular_formula.rstrip('.+[]') + ']+'
                                    molecular_formula_anion = '[' + molecular_formula.rstrip('.-+[]') + ']-'
                                    
                                    # remove any . after the masses
                                    calculated_mass_from_si = str_split[line_structure['calculated mass from si']]
                                    calculated_mass_from_si = calculated_mass_from_si.rstrip('.')
                                    calculated_mass_from_si = float("{:.4f}".format(float(calculated_mass_from_si)))

                                    found_mass_from_si = str_split[line_structure['found mass from si']]
                                    found_mass_from_si = found_mass_from_si.rstrip('.')
                                    found_mass_from_si = float("{:.4f}".format(float(found_mass_from_si)))
                                    
                                    # calculate the masses based on the molecular formula
                                    calculated_mass_from_neutral = float("{:.4f}".format(Formula(molecular_formula).monoisotopic_mass))
                                    calculated_mass_from_cation = float("{:.4f}".format(Formula(molecular_formula_cation).monoisotopic_mass))
                                    calculated_mass_from_anion = float("{:.4f}".format(Formula(molecular_formula_anion).monoisotopic_mass))

                                    # determine measuring mode here:
                                    initial_comment = ""
                                    # molecular_ion_type = ""
                                    if calculated_mass_from_neutral == calculated_mass_from_si: # neutral mode
                                        # molecular_ion_type = "neutral"
                                        measuring_mode = "neutral"
                                    elif calculated_mass_from_cation == calculated_mass_from_si: # cation mode
                                        # molecular_ion_type = "cation"
                                        measuring_mode = "cation"
                                    elif calculated_mass_from_anion == calculated_mass_from_si: # anion mode
                                        # molecular_ion_type = "anion"
                                        measuring_mode = "anion"
                                    else:
                                        # molecular_ion_type = "unknown"
                                        measuring_mode = "cation"


                                    extracted_data.loc[len(extracted_data)] = [page_num+1, found_string, initial_comment, molecular_ion_type, molecular_formula, molecular_formula_cation, molecular_formula_anion, 
                                                                            measuring_mode, sodium, calculated_mass_from_si, found_mass_from_si, 
                                                                            calculated_mass_from_neutral, calculated_mass_from_cation, calculated_mass_from_anion, 
                                                                            calculateError(found_mass_from_si, calculated_mass_from_si), 
                                                                            calculateError(found_mass_from_si, calculated_mass_from_neutral), calculateError(found_mass_from_si, calculated_mass_from_cation), calculateError(found_mass_from_si, calculated_mass_from_anion)]
                                curr_index += 100
                        except (ValueError, FormulaError) as e:
                            match e:
                                case ValueError():
                                    print("val error")
                                    extracted_data.loc[len(extracted_data)] = [page_num+1, str_split, 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A']
                                case FormulaError():
                                    print("formula error")
                                    # print(str_split)
                                    calculated_mass_from_si = str_split[line_structure['calculated mass from si']]
                                    calculated_mass_from_si = calculated_mass_from_si.rstrip('.')
                                    calculated_mass_from_si = float("{:.4f}".format(float(calculated_mass_from_si)))

                                    found_mass_from_si = str_split[line_structure['found mass from si']]
                                    found_mass_from_si = found_mass_from_si.rstrip('.')
                                    found_mass_from_si = float("{:.4f}".format(float(found_mass_from_si)))

                                    extracted_data.loc[len(extracted_data)] = [page_num+1, str_split, -1, molecular_ion_type, str_split[line_structure['molecular formula']], -1, -1, 
                                        -1, -1, calculated_mass_from_si, found_mass_from_si, 
                                        -1, -1, -1, 
                                        calculateError(found_mass_from_si, calculated_mass_from_si), 
                                        -1, -1, -1]

                            
                            # print(e)
                            curr_index += 100
                            continue
                    

            except ValueError as e: # on this page, there is no hrms data
                continue
        
        print(len(extracted_data))
        # for every entry in the pandas df
        for index, row in extracted_data.iterrows(): # page number = row['page number']
            # print(row['molecular formula neutral'])
            comment = row['initial comment']
            errlvl = ""

            total_examples = total_examples + 1
            

            if isinstance(row['measuring mode'], int) and row['measuring mode'] == -1:
                comment = f"Invalid molecular formula. Check for capitalizations, notations (i.e. 0's mistaken for O's), or any other typographical errors. "
                errlvl = "B"
                incorrect_examples[1] = incorrect_examples[1] + 1
            elif isinstance(row['measuring mode'], str) and row['measuring mode'] == 'N/A':
                comment = f"The SI provided could not be processed due to unexpected inputs. "
                invalid_inputs += 1
            elif row['mass error from ' + row['measuring mode']] >= threshold or row['mass error from si'] >= threshold or row['measuring mode'] == 'neutral' or row['molecular ion type'] == 'unknown':
                try:
                    if row['mass error from ' + row['measuring mode']] > 100 or row['mass error from si'] > 100:
                        composition = Formula(row['molecular formula neutral']).composition()
                        for element, content in composition.items():
                            if content.count >= 100:
                                comment = "Potential invalid molecular formula. Check for capitalizations, notations (i.e. 0's mistaken for O's), or any other typographical errors. "
                                errlvl = "B"
                            if 'I' in element or 'l' in element:
                                comment = "Potential invalid molecular formula. Check for capitalizations, notations (i.e. 0's mistaken for O's), or any other typographical errors. "
                                errlvl = "B"
                    if row['measuring mode'] == 'neutral':
                        if not neutral: 
                            comment += "The reported mass was calculated not taking into account the mass of the electron. "
                            errlvl = "F"
                        elif neutral:
                            comment += "Above selected threshold. "
                            errlvl = "G"
                    
                    molecular_weight_neutral = float("{:.4f}".format(Formula(row['molecular formula neutral']).mass))
                    molecular_weight_ion = float("{:.4f}".format(Formula(row['molecular formula ' + row['measuring mode']]).mass))

                    if row['calculated mass from si'] == molecular_weight_neutral:
                        comment += f"The molecular weight ({molecular_weight_neutral}) was calculated, not the accurate mass ({row['calculated mass from ' + row['measuring mode']]}). "
                        errlvl = "C"
                    if row['calculated mass from si'] == molecular_weight_ion:
                        comment += f"The molecular weight ({molecular_weight_ion}) was calculated, not the accurate mass ({row['calculated mass from ' + row['measuring mode']]}). "
                        errlvl = "C"
                    
                    # check H, Na, K for neutral-actualmass, neutral, mw_neutral-actualmass, mw_neutral
                    ions_to_check = {
                        'H': 1,
                        'Na': 23,
                        'K': 39
                    }
                    for element in ions_to_check:
                        assumed_mass = ions_to_check[element]

                        neutral_assumed = float("{:.4f}".format(row['calculated mass from neutral'] + assumed_mass))
                        neutral_actual_assumed = float("{:.4f}".format(row['calculated mass from neutral'] - Formula(element).monoisotopic_mass + assumed_mass))
                        mw_assumed = float("{:.4f}".format(molecular_weight_neutral + assumed_mass))
                        mw_actual_assumed = float("{:.4f}".format(molecular_weight_neutral - Formula(element).mass + assumed_mass))

                        composition = compositionToDict(Formula(row['molecular formula neutral']).composition())

                        if neutral_assumed == row['calculated mass from si']:
                            composition[element]['count'] = composition[element]['count'] + 1
                            molecular_formula_addition = compositionToFormula(composition, 'neutral')
                            comment += f"It appears that the accurate mass was generated by calculating the accurate mass for the neutral molecule {row['molecular formula neutral']} ({row['calculated mass from neutral']}) and adding +{assumed_mass}.0000 => {neutral_assumed}. "
                            errlvl = "E"
                        elif element in composition and neutral_actual_assumed == row['calculated mass from si']:
                            composition[element]['count'] = composition[element]['count'] - 1
                            molecular_formula_removed = compositionToFormula(composition, 'neutral')
                            comment += f"It appears that the accurate mass was generated by calculating the accurate mass for the neutral molecule {molecular_formula_removed} ({Formula(molecular_formula_removed).monoisotopic_mass:.4f}) and adding +{assumed_mass:.4f} => {row['molecular formula neutral']} ({neutral_actual_assumed}). "
                            errlvl = "E"
                        elif mw_assumed == row['calculated mass from si']:
                            composition[element]['count'] = composition[element]['count'] + 1
                            molecular_formula_addition = compositionToFormula(composition, 'neutral')
                            comment += f"It appears that the accurate mass was generated by calculating the molecular weight for the neutral molecule {row['molecular formula neutral']} ({molecular_weight_neutral}) and adding +{assumed_mass}.0000 => {mw_assumed}. "
                            errlvl = "C"
                        elif element in composition and mw_actual_assumed == row['calculated mass from si']:
                            composition[element]['count'] = composition[element]['count'] - 1
                            molecular_formula_removed = compositionToFormula(composition, 'neutral')
                            comment += f"It appears that the accurate mass was generated by calculating the molecular weight for the neutral molecule {molecular_formula_removed} ({Formula(molecular_formula_removed).mass:.4f}) and adding +{assumed_mass:.4f} => {row['molecular formula neutral']} ({mw_actual_assumed:.4f})."
                            errlvl="C"
                    # calc_found = row['calculated mass from si'] - row['found mass from si']
                    # if errlvl == "" and (abs(row['calculated mass from si'] - row['calculated mass from ' + row['measuring mode']]) > calc_found or abs(row['found mass from si'] - row['calculated mass from ' + row['measuring mode']]) > calc_found):
                    if errlvl == "" and row['mass error from si'] < row['mass error from ' + row['measuring mode']] - abs(row['mass error from si'] < row['mass error from ' + row['measuring mode']]) > 1:
                        comment += "Found mass matches erroneous formula, recheck data and/or look for extraneous or missing atom(s) in the reported molecular formula. "
                        errlvl = "A"

                    
                    if errlvl == "":
                        # check swapped
                        if checkSwappedIntegers(row['calculated mass from si'], row['found mass from si']):
                            comment += "The calculated mass and the measured mass appear to be transposed by two digits. "
                            errlvl = "D"
                        
                        # check diff. integers. does this need to be last?
                        elif checkDifferentIntegers(row['calculated mass from si'], row['found mass from si']):
                            comment += "The reported and measured accurate masses differ in their integers. "
                            errlvl = "D"
                        # check typo between si calculated/found
                        elif checkTypo(row['calculated mass from si'], row['found mass from si'], row['calculated mass from ' + row['measuring mode']]):
                            comment += "The calculated mass might contain a typo. "
                            errlvl = "D"
                    
                    # determine whether C
                    
                    if errlvl == "" and threshold < 5 and row['mass error from si'] <= 5 and row['mass error from ' + row['measuring mode']] <= 5:
                        comment += "Above selected threshold. "
                        errlvl = "G"
                    elif errlvl == "" and (row['mass error from si'] > threshold or row['mass error from ' + row['measuring mode']]):
                        comment += "Found mass matches erroneous formula, recheck data and/or look for extraneous or missing atom(s) in the reported molecular formula. "
                        errlvl = "A"
                    

                    if errlvl == "A":
                        incorrect_examples[0] = incorrect_examples[0] + 1
                    elif errlvl == "F":
                        incorrect_examples[2] = incorrect_examples[2] + 1
                    elif errlvl == "G":
                        incorrect_examples[3] = incorrect_examples[3] + 1
                    elif errlvl != "":
                        incorrect_examples[1] = incorrect_examples[1] + 1

                    # # check adding ions
                    # ions_to_check = {
                    #     'H': 1,
                    #     'NH4': 18,
                    #     'Na': 23,
                    #     'K': 39,
                    #     'Cl': 35,
                    #     'CH3COO': 59
                    # }
                    # for ion in ions_to_check:
                    #     assumed_mass = ions_to_check[ion]
                    #     check_ion_result = checkIon(row['molecular formula neutral'], row['calculated mass from si'], row['calculated mass from neutral'], ion, assumed_mass)
                    #     print(check_ion_result)
                    #     if check_ion_result == 1:
                    #         comment += f"It appears that the high resolution mass was generated by calculating the accurate mass for the neutral molecule {row['molecular formula neutral']} ({row['calculated mass from neutral']}) and adding +{assumed_mass}.0000 => ({row['calculated mass from neutral']+assumed_mass})"
                    #     elif check_ion_result == 2:
                    #         molecular_formula_neutral_composition = compositionToDict(Formula(row['molecular formula neutral']).composition())
                    #         ion_composition = compositionToDict(Formula(ion).composition())
                    #         for element in ion_composition:
                    #             print(element)
                    #             molecular_formula_neutral_composition[element]['count'] = molecular_formula_neutral_composition[element]['count'] - ion_composition[element]['count']
                    #         # molecular_formula_neutral_composition[ion]['count'] = molecular_formula_neutral_composition[ion]['count'] - 1
                    #         molecular_formula_neutral_reduced = compositionToFormula(molecular_formula_neutral_composition, "neutral")
                    #         calculated_mass_from_neutral_reduced = float("{:.4f}".format(Formula(molecular_formula_neutral_reduced).monoisotopic_mass))
                    #         comment += f"It appears that the high resolution mass was generated by calculating the accurate mass for the neutral molecule {molecular_formula_neutral_reduced} ({calculated_mass_from_neutral_reduced}) and adding +{assumed_mass}.0000 => ({calculated_mass_from_neutral_reduced+assumed_mass})"
                    # print("finish ions")
                    # composition = compositionToDict(Formula(row['molecular formula ' + row['measuring mode']]).composition())
                    # # print(composition)
                    # for element in composition:
                    #     if element != 'e-':
                    #         # check addition/remove
                    #         check_addition = checkAddition(row['molecular formula ' + row['measuring mode']], row['calculated mass from si'], element, 10)
                    #         check_remove = checkRemove(row['molecular formula ' + row['measuring mode']], row['calculated mass from si'], element, 10)
                    #         if check_addition != 0:
                    #             composition[element]['count'] = composition[element]['count'] + check_addition
                    #             new_molecular_formula = compositionToFormula(composition, row['measuring mode'])                            
                    #             comment += f"Adding {check_addition} {element}-atom(s), the molecular formula fits the mass reported in the SI: {new_molecular_formula} ({Formula(new_molecular_formula).monoisotopic_mass:.4f})"
                    #         if check_remove != 0:
                    #             composition[element]['count'] = composition[element]['count'] - check_remove
                    #             new_molecular_formula = compositionToFormula(composition, row['measuring mode'])                            
                    #             comment += f"Removing {check_remove} {element}-atom(s), the molecular formula fits the mass reported in the SI: {new_molecular_formula} ({Formula(new_molecular_formula).monoisotopic_mass:.4f})"

                    # composition = compositionToDict(Formula(row['molecular formula ' + row['measuring mode']]).composition())
                    # # check addition for elements not in the molecular formula
                    # elements_to_check_addition = ['D', 'H', 'Li', 'B', 'C', 'O', 'F', 'Na', 'Si', 'P', 'S', 'Cl', 'K', 'Br', 'I']
                    # for element in elements_to_check_addition: 
                    #     if element not in composition:
                    #         check_addition = checkAddition(row['molecular formula ' + row['measuring mode']], row['calculated mass from si'], element, 5)
                    #         if check_addition != 0:
                    #             composition[element] = {"count": check_addition}
                    #             new_molecular_formula = compositionToFormula(composition, row['measuring mode'])     
                    #             comment += f"Adding {check_addition} {element}-atom(s), the molecular formula fits the mass reported in the SI: {new_molecular_formula} ({Formula(new_molecular_formula).monoisotopic_mass:.4f})"                        

                    # composition = compositionToDict(Formula(row['molecular formula ' + row['measuring mode']]).composition())
                    # # check addition/remove for element groups
                    # element_groups_to_check = ['CH2', 'CH3', 'CH4', 'OH', 'H2O', 'H3O', 'NH', 'NH2', 'NH3', 'NH4']
                    # for element in element_groups_to_check:
                    #     # check addition/remove
                    #     check_addition = checkAddition(row['molecular formula ' + row['measuring mode']], row['calculated mass from si'], element, 5)
                    #     check_remove = checkRemove(row['molecular formula ' + row['measuring mode']], row['calculated mass from si'], element, 5)
                    #     if check_addition != 0:
                    #         element_composition = compositionToDict(Formula(element).composition())
                    #         for curr_element in element_composition:
                    #             if curr_element not in composition:
                    #                 composition[curr_element] = {"count": check_addition * element_composition[curr_element]['count']}
                    #             else:
                    #                 composition[curr_element]['count'] += check_addition * element_composition[curr_element]['count']
                    #         new_molecular_formula = compositionToFormula(composition, row['measuring mode'])                            
                    #         comment += f"Adding {check_addition} {element}-atom(s), the molecular formula fits the mass reported in the SI: {new_molecular_formula} ({Formula(new_molecular_formula).monoisotopic_mass:.4f})"
                    #     if check_remove != 0:
                    #         element_composition = compositionToDict(Formula(element).composition())
                    #         for curr_element in element_composition:
                    #             composition[curr_element]['count'] += -check_remove * element_composition[curr_element]['count']
                    #         new_molecular_formula = compositionToFormula(composition, row['measuring mode'])      
                    #         comment += f"Removing {check_remove} {element}-atom(s), the molecular formula fits the mass reported in the SI: {new_molecular_formula} ({Formula(new_molecular_formula).monoisotopic_mass:.4f})"

                    # composition = compositionToDict(Formula(row['molecular formula ' + row['measuring mode']]).composition())

                    # # check replace H+ -> Na+, Na+ -> H+, Na+ -> K+, K+ -> Na+
                    # print("start of replace")
                    # replace_elements_dictionary = {
                    #     'H' : ['Na', 'K', 'NH4', 'Cl', 'CH3COO'],
                    #     'Na': ['H', 'K', 'NH4'],
                    #     'K': ['H', 'Na', 'NH4'],
                    #     'Cl': ['H', 'CH3COO'],
                    #     'CH3COO': ['H', 'Cl']
                    # }
                    # for element in replace_elements_dictionary:
                    #     if element in composition:
                    #         value = replace_elements_dictionary[element]
                    #         for curr_value in value:
                    #             check_replace_neutral = checkReplace(row['molecular formula ' + row['measuring mode']], row['calculated mass from si'], element, curr_value) 
                    #             if check_replace_neutral:
                    #                 composition[element]['count'] = composition[element]['count'] - 1
                    #                 if curr_value not in composition:
                    #                     composition[curr_value] = {"count": 1}
                    #                 else:
                    #                     composition[curr_value]['count'] = composition[curr_value]['count'] + 1
                    #                 new_molecular_formula = compositionToFormula(composition, row['measuring mode'])
                    #                 comment += f"Replacing {element} with {curr_value}, the molecular formula fits the mass reported in the SI: {new_molecular_formula} ({Formula(new_molecular_formula).monoisotopic_mass:.4f})"
                                        
                except ValueError as e:
                    print(e)
                    comment = f"On page {row['page number']}, found invalid line, '{' '.join(row['extracted text'])}'."
            
            # determine error level 
            if row['measuring mode'] != -1 and row['measuring mode'] != 'N/A':
                file_request.append({
                    "molform": row['molecular formula neutral'],
                    "pg": row['page number'],
                    "iontype": row['molecular ion type'],
                    "errlvl": errlvl,
                    "errms": row['mass error from si'],
                    "errcalc": row['mass error from ' + row['measuring mode']],
                    "sicalc": row['calculated mass from si'],
                    "sifound": row['found mass from si'],
                    "recalc": row['calculated mass from ' + row['measuring mode']],
                    "com": comment
                })
            elif row['measuring mode'] == 'N/A':
                file_request.append({
                    "molform": 'N/A',
                    "pg": row['page number'],
                    "iontype": 'N/A',
                    "errlvl": 'N/A',
                    "errms": 'N/A',
                    "errcalc": 'N/A',
                    "sicalc": 'N/A',
                    "sifound": 'N/A',
                    "recalc": 'N/A',
                    "com": comment
                })
            else:
                # print(errlvl)
                # determine if invalid formula is a POTENTIAL or KNOWN error ?

                if row['mass error from si'] != -1:
                    err_ms = row['mass error from si']
                else:
                    err_ms = "N/A"

                file_request.append({
                    "molform": row['molecular formula neutral'],
                    "pg": row['page number'],
                    "iontype": row['molecular ion type'],
                    "errlvl": errlvl,
                    "errms": err_ms,
                    "errcalc": "N/A",
                    "sicalc": row['calculated mass from si'],
                    "sifound": row['found mass from si'],
                    "recalc": "N/A",
                    "com": comment
                })

    except Exception as e:
        print("Error: ", e)
        continue

    # convert the data into json format + write the json data into a file
    # the file will be read by the server later to generate the table

    file_json = {
        "Title": title,
        "filepath": filepath,
        "Date": date.today().isoformat(),
        "threshold": threshold,
        "total": total_examples,
        "aerrors": incorrect_examples[0],
        "bgerrors": incorrect_examples[1],
        "herrors": incorrect_examples[2],
        "ierrors": incorrect_examples[3],
        "invalidInputs": invalid_inputs,
        "file_request": file_request
    }

    json_data = json.dumps(file_json, indent=3)
    outputJSONFile.write(json_data)
    outputJSONFile.write("qwqwqw\n")

outputJSONFile.close()


