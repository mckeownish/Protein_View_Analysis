# Library Imports

import streamlit as st
import CarbonaraDataTools as CDT
import biobox
import glob

from stmol import showmol,obj_upload
import py3Dmol

import numpy as np
import subprocess


tab1, tab2, tab3 = st.tabs(["PDB_viewer", "Carbonara_run", "Carbonara_results"])


# st.title('CARBONARA Script Builder')


with tab1:
    st.title('CARBONARA Visualisation')

# -*-*-*- PDB viewer -*-*-*-

    uploaded_file = st.sidebar.file_uploader("Upload PDB file")

    if uploaded_file is not None:
        obj = obj_upload(uploaded_file)
        showmol(obj,width=800)

    # -------------------------------- 


    st.write('To implement drag + drop - issues with byte vs string objects in DSSP prediction..')
    # pdb_file_upload = st.file_uploader('Upload a PDB file!', type='pdb')


    


    # if pdb_file_upload is not None:
    # #     pdb_file = pdb_file.decode('utf-8')
    #     with open("temp_file.pdb", 'w') as temp:
    #         temp.write(pdb_file_upload.getvalue().decode('utf-8'))


    # st.write(ss)


with tab2:

    st.title('CARBONARA Runner')

    # -*-*-*- Writing files in correct format + expected locations -*-*-*-

    pdb_file = 'Example/Lysozyme/Lysozyme.pdb'  
    SAXS_file = 'Example/Lysozyme/LysozymeSaxs.dat'

    working_path = CDT.setup_working_directory()

    M = CDT.pdb_2_biobox(pdb_file)

    # Extract coordinates + primary sequence
    coords = CDT.extract_CA_coordinates(M)
    # st.write(coords)

    sequence = CDT.extract_sequence(M)

    secondary = CDT.DSSP_structure_extractor(pdb_file)

    CDT.write_fingerprint_file(1, sequence, secondary, working_path)
    CDT.write_coordinates_file(coords, working_path)
    CDT.write_mixture_file(working_path)

    CDT.write_saxs(SAXS_file, working_path)


    # -*-*-*- Writing expected format file -*-*-*-

    run_auto = True

    @st.cache_data
    def auto_varying_function(run_auto):
        if run_auto:
            return CDT.find_non_varying_linkers()


    allowed_linker, linker_indices = auto_varying_function(run_auto)


    selected_linkers_lst = st.multiselect('Select allowed varying linkers',
                                     linker_indices,
                                     allowed_linker,
                                     label_visibility="collapsed")

    # st.write(selected_linkers_lst)

    q_exp_min, q_exp_max, q_Q1, q_Q3 = CDT.get_minmax_q(SAXS_file)


    min_q, max_q = list( st.slider( 'Select a range of q values', q_exp_min, q_exp_max, (q_exp_min, q_Q3), step=0.01 ))
    fig = CDT.SAXS_selection_plotter(SAXS_file, min_q, max_q)

    st.plotly_chart(fig)

    col1, col2 = st.columns(2)

    fit_n_times = col1.number_input('Choose number of Carbonara repeats', min_value=1, value=5)
    max_fit_steps = col2.number_input('Choose naximum number of fit steps in each repeat', min_value=10, value=5000, step=100)
    
    @st.cache_data
    def run_commandline_carbonara(execute):
        if execute:

            CDT.write_varysections_file(selected_linkers_lst, working_path)
            CDT.write_sh_file(working_path=working_path, fit_n_times=fit_n_times, min_q=min_q, max_q=max_q, max_fit_steps=max_fit_steps)
            st.write('sh file written!')
            result = subprocess.run(['bash', 'RunMe.sh'], capture_output=True, text=True)

    
    execute = st.button('Run Carbonara')

    if execute:
        run_commandline_carbonara(execute)
        
        st.write("WOOO, you've run Carbonara entirely code-freeee babyyy!")
        # st.write(result.stdout)


with tab3:

    st.title('CARBONARA Analysis')

    run_selection = st.selectbox('Select run number', np.arange(fit_n_times)+1)

    scattering_lst = CDT.sort_by_creation(glob.glob('Fitting/fitdata/fitmolecule'+str(run_selection)+'SubstepScatter*'))
    substruct_lst = CDT.sort_by_creation(glob.glob('Fitting/fitdata/fitmolecule'+str(run_selection)+'Substep_*'))

    if len(scattering_lst) != len(substruct_lst): 
        st.write('uh oh, your scattering and structure output lists are different lengths :(')


    sub_selection = st.selectbox('Select substep', np.arange(len(substruct_lst))+1)