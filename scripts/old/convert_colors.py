 
import matplotlib.colors

colors = {
"#228b22":"forestgreen",
"#00008b":"darkblue",
"#b03060":"maroon3",
"#ff4500":"orangered",
"#ffff00":"yellow",
"#deb887":"burlywood",
"#00ff00":"lime",
"#00ffff":"aqua",
"#ff00ff":"fuchsia",
"#6495ed":"cornflower",
} # from https://mokole.com/palette.html

for c in colors:
	rgb = matplotlib.colors.to_rgb(c)
	rgb = ' '.join("%01.3f"%x for x in rgb)
	print(c,"\t",rgb,"\t",colors[c],)