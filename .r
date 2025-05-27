ft=7;
        w=14;h=7;
        o="g.png";
opts="14";
        library(ComplexHeatmap);
        tt=read.table("stdin",header=T,sep="\t");
        m=as.matrix(tt[,-1]);
        row.names(m)=tt[,1];
        h=Heatmap(m,column_names_gp=gpar(fontsize=ft), row_names_gp = gpar(fontsize = ft) );
        png(o,width=w,height=h)
