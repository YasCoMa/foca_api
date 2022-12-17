
import pandas as pd
import pickle

#Sequences by states 2 month
@app.route('/get_statesbr')
def section14_get_statesbr ():
    stbr = pd.read_pickle ("/home/tavinbio/public_html/foca_backend/data/Numbers_sequences_Brazil_states.pkl")
    jstbr = stbr.to_json (orient='records')
    return json.dumps (jstbr)
