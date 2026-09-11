#!/usr/bin/env python3
from pathlib import Path
import re
HERE=Path(__file__).resolve().parent

def module_global(filename, global_name, exports):
    source=(HERE.parent/'web'/filename).read_text(encoding='utf-8')
    source=re.sub(r'^export ', '', source, flags=re.MULTILINE)
    return '<script>\nconst '+global_name+'=(()=>{\n'+source+'\nreturn '+exports+';})();\n</script>'

def main():
    modules={
      'model':module_global('model.mjs','Wilberforce','{validate,normalModes,analyticState,rk4Step,energyOf,advanceClock}'),
      'presets':module_global('presets.mjs','WilberforcePresets','PRESETS')}
    for name in ('index','index_lite'):
        text=(HERE/(name+'.template.html')).read_text(encoding='utf-8')
        def inline(match):
            key=match.group(1)
            if key in modules:return modules[key]
            return '<script>\n'+(HERE/key).read_text(encoding='utf-8')+'\n</script>'
        text=re.sub(r'<!--INLINE:([^>]+)-->',inline,text)
        if name=='index':
            notice=(HERE/'vendor/LICENSE-three.txt').read_text(encoding='utf-8')
            text=text.replace('<head>','<head>\n<!--\n'+notice+'\n-->')
        (HERE/(name+'.html')).write_text(text,encoding='utf-8')
        print('wrote',name+'.html',len(text),'characters')
if __name__=='__main__':main()
