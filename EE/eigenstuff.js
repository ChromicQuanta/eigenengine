function findRealEigenVector(matr,eigenValue,tol){
var mat = cloneMat(matr)
mat = idmin(eigenValue,mat)
for(let n=1;n<mat.length;n++){

    for(let m=0;m<n;m++){
        var comp = cmult(dotp(mat[n],mat[m]),crpow(vmag(mat[m]),-2))
        //console.log(comp.r)
        if(cmult(comp,cconj(comp)).r>(1/tol)) continue
        
mat[n]=subRows(mat[n],cmultRow(comp,mat[m]))
}
    if(isZero(mat[n],tol)){
        break
    }
}

var rand = []
for(let n=0;n<mat.length;n++){
    rand[n] = z(Math.random()*2-1,Math.random()*2-1)
}

for(let n=0;n<mat.length;n++){
    //console.log(rand)
        if(isZero(mat[n],tol)){
        continue
    }
rand=subRows(rand,cmultRow(cmult(dotp(rand,mat[n]),crpow(vmag(mat[n]),-2)),mat[n]))

}
return rand
}
function measure(state){
    var left = 1
    for(let n=0;n<state.length;n++){

var prob = Math.hypot(state[n][0].r,state[n][0].i)**2
var ret = Math.random()<=(prob/left)     
        left-=prob

        if(ret) return n
    }
}

function det(matr){
   var mat = cloneMat(matr)
    for(let n=1;n<mat.length;n++){
        for(let m=0;m<n;m++){
            if(mat[m][m].r==0 && mat[m][m].i==0){
                for(let o=0;o<mat.length;o++){
                    if(mat[o][m].r!=0 || mat[o][m].i!=0){
                        mat[m]=addRows(mat[o],mat[m])
                        break
                }
                    
                if(o==mat.length-1) return z(0,0)
            }
            }
            mat[n]=addRows(mat[n],cmultRow(cmult(csub(z(0,0),mat[n][m]),cinv(mat[m][m])),mat[m]))
        }
    }
    prod = z(1,0)
    for(let n=0;n<mat.length;n++){
        prod = cmult(mat[n][n],prod)
    }
    return prod
}








function charpoly(matr){
mat = cloneMat(matr)

var vec = []
        for(n=0;n<mat.length+1;n++){
        vec[n] = det(idmin(1+(10*n/(matr.length+1)),mat))
    }

    var mat = powermat(vec.length)
    for(let n=1;n<mat.length;n++){
    for(let m=0;m<n;m++){
        var c = cmult(mat[n][m],cinv(mat[m][m]))
        vec[n] = csub(vec[n],cmult(c,vec[m]))
        mat[n] = subRows(mat[n],cmultRow(c,mat[m]))
    }
    }

    for(let n=mat.length-2;n>=0;n--){
    for(let m=mat.length-1;m>n;m--){
        var c = cmult(mat[n][m],cinv(mat[m][m]))
        vec[n] = csub(vec[n],cmult(c,vec[m]))
        mat[n] = subRows(mat[n],cmultRow(c,mat[m]))
    }
    }
    for(let n=0;n<vec.length;n++){
        vec[n] = cmult(cinv(mat[n][n]),vec[n])
    }
    return vec
}


function eigenSolveHamiltonian(mat,delta,gap,limit,tol){
 roots = findRealRoots(charpoly(mat),delta,gap,limit)
var eigenVectors = []
for(let n=0;n<roots.length;n++){
var eigen = findRealEigenVector(mat,roots[n].r,tol)
if(!isZero(eigen,tol)){
eigenVectors.push(normalise(eigen))
}
}
return eigenVectors
}





function realRootFinder(poly,start,delta,limit){
    del = delta
    var pt = z(start,0)
    var der = derPoly(poly)

    var lastSign = 1
    for(let n=0;n<limit;n++){
    var p = evalPoly(poly,pt)
    var d = evalPoly(der,pt)

    pt.r-=Math.sign(d.r)*cap1(p.r)*del
    p = evalPoly(poly,pt)
    if(n && lastSign!=Math.sign(p.r)) del/=10
    lastSign=Math.sign(p.r)
    //pt.i-=d.i*cap1(p.i)*delta
    }
    return [pt,p]
}
function findRealRoots(poly,delta,gap,limit){
    var polys=[]
    var rt=[[],[]]
    polys[0] = poly
    for(let n=0;n<poly.length-2;n++){
        polys[n+1]=derPoly(polys[n])
    }

    rt[0]=[ cmult(csub(z(0,0),polys[polys.length-1][0]),cinv(polys[polys.length-1][1])) ]
    //return [polys[polys.length-1],roots]
    for(let n=polys.length-2;n>=0;n--){
        rt[1]=[]
       // console.log(rt[0])
        for(let m=0;m<rt[0].length;m++){
            if(m==0 || 0){
            rt[1].push(realRootFinder(polys[n],rt[0][m].r-gap,delta,limit)[0])
            rt[1].push(realRootFinder(polys[n],rt[0][m].r+gap,delta,limit)[0])
            }else{
            rt[1].push(realRootFinder(polys[n],rt[0][m].r+gap,delta,limit)[0])
            }
        }
        
        rt.reverse()
    }
    return rt[0]
}