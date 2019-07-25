import requests
import json
import pandas as pd

from chimera import df_pfam


if __name__ == '__main__':

    pfam_names = df_pfam['pfam_id'].unique()

    headers = {'Content-type': 'application/json', 'Accept': 'application/json'}
    data = {
        "hmmdb": "pfam",
        "cut_ga": True,
        "seq": "MEGDAVEAIVEESETFIKGKERKTYQRRREGGQEEDACHLPQNQTDGGEVVQDVNSSVQMVMMEQLDPTLLQMKTEVMEGTVAPEAEAAVDDTQIITLQVVNMEEQPINIGELQLVQVPVPVTVPVATTSVEELQGAYENEVSKEGLAESEPMICHTLPLPEGFQVVKVGANGEVETLEQGELPPQEDPSWQKDPDYQPPAKKTKKTKKSKLRYTEEGKDVDVSVYDFEEEQQEGLLSEVNAEKVVGNMKPPKPTKIKKKGVKKTFQCELCSYTCPRRSNLDRHMKSHTDERPHKCHLCGRAFRTVTLLRNHLNTHTGTRPHKCPDCDMAFVTSGELVRHRRYKHTHEKPFKCSMCDYASVEVSKLKRHIRSHTGERPFQCSLCSYASRDTYKLKRHMRTHSGEKPYECYICHARFTQSGTMKMHILQKHTENVAKFHCPHCDTVIARKSDLGVHLRKQHSYIEQGKKCRYCDAVFHERYALIQHQKSHKNEKRFKCDQCDYACRQERHMIMHKRTHTGEKPYACSHCDKTFRQKQLLDMHFKRYHDPNFVPAAFVCSKCGKTFTRRNTMARHADNCAGPDGVEGENGGETKKSKRGRKRKMRSKKEDSSDSENAEPDLDDNEDEEEPAVEIEPEPEPQPVTPAPPPAKKRRGRPPGRTNQPKQNQPTAIIQVEDQNTGAIENIIVEVKKEPDAEPAEGEEEEAQPAATDAPNGDLTPEMILSMMDR"
    }
    r = requests.post('https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan', data=json.dumps(data), headers=headers).json()
    hits = r['results']['hits']

    l = []
    for hit in hits:
        for d in hit['domains']:
            pfam_name = d['alihmmacc'][:7] + '_' + d['alihmmname']
            l.append({
                'pfam_domain': pfam_name,
                'target_start': int(d['alisqfrom']),
                'target_end': int(d['alisqto']),
                'hmm_start': int(d['alihmmfrom']),
                'hmm_end': int(d['alihmmto']),
                'domain_length': int(d['aliM']),
                'bit_score': float(d['bitscore']),
                'reported': bool(d['is_reported']),
                'e_value': float(d['ievalue']),
                'aliseq': d['aliaseq'],
                'interacdome': pfam_name in pfam_names
            })

    df = pd.DataFrame(l)

    df = df[(df['bit_score']>0) & (df['hmm_start']==1) & (df['hmm_end']==df['domain_length'])]
    df = df.sort_values('target_start')
    print(df)
