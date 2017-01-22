
for shake in 15 25 35 45 55 65; do 
	for linker in 100 120 140 160 180 200; do
		for i in {1..10}; do
			python process.py l${linker}_s${shake}_${i} $linker $shake
		done
	done
done
