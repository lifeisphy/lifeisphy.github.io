import sys,pywintypes

rep={
'https://cdn.jsdelivr.net/npm/katex@0.15.3/dist/katex.min.css':
'https://gkone-pku-phy.github.io/assets/dependencies/katex/katex.min.css',

'https://cdn.jsdelivr.net/npm/reveal.js@4.1.0/dist/theme/white.css': 

'https://gkone-pku-phy.github.io/assets/dependencies/reveal/css/theme/white.css',
      
'https://cdn.jsdelivr.net/npm/reveal.js@4.1.0/dist/reveal.js':
'https://gkone-pku-phy.github.io/assets/dependencies/reveal/js/reveal.js'
}
file = sys.argv[1]
with open(file,'r') as f:
    s:str=f.read()
for k,v in rep.items():
    s=s.replace(k,v)
with open(file,'w') as f:
    f.write(s)
