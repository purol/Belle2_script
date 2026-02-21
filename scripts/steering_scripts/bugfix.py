with open('./gbasf2_Kumoi_BUGFIX3_CWR.py', 'rb+') as f:
    content = f.read()
    f.seek(0)
    f.write(content.replace(b'\r', b''))
    f.truncate()
