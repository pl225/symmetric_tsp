dir=instancias/
for arquivo in `ls $dir`
    do
        #echo "INSTANCIA"
        #echo "$arquivo"
        echo "$(./exe "${dir}$arquivo")"
done