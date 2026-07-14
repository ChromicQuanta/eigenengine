function logmat(mat,pad,round){
    var n
    var m
var ans = ""
    for(n=0;n<mat.length;n++){
    for(m=0;m<mat[0].length;m++){
        ans+="("+cround(mat[n][m],round).string.padStart(pad," ")+")"
    }
ans+="\n"
    }
    console.log(ans)
}
function powermat(dim){
    var ans = []
    
    for(let n=0;n<dim;n++){
        ans[n]=[]
    for(let m=0;m<dim;m++){
            ans[n][m]=z((1+10*(n/dim))**m,0)
    } 
    }
    return ans
}
function construct1DHamiltonian(dim,Vfunc,x1,x2){
    ans = []
    for(var n=0;n<dim;n++){
        ans[n]=[]
        for(var m=0;m<dim;m++){
        if(n==m){
            ans[n][m]=z((-2+Vfunc(x1+(x2-x1)*n/(dim-1))),0)
            continue
        }
        if(n==(m-1) || (n-1)==m){
            ans[n][m]=z(1,0)
            continue
        }
        ans[n][m]=z(0,0)
    }
    }
    //logmat(ans,5,3)
return ans
}