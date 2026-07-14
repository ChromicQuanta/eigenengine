    function ComplexNumber(real,imag){
    this.r = real
    this.i = imag
    Object.defineProperty(this,"string",{
        get:(x)=>{var x = this;
                  function i1(i,m){
                      if(i==1) return ""
                      if(i==-1 && !m) return "-"
                      if(i==-1 && m) return ""
                      return (m?-1:1)*i+""
                  }
            if(x.i==0)return x.r+""
            if(x.r==0)return i1(x.i)+"i"
            if(x.i>0)return x.r + " + " + i1(x.i)+"i"
            if(x.i<0)return x.r + " - " + i1(x.i,1)+"i"
            
            
        }
    })
}
function isZero(vec,tol){
    for(let n=0;n<vec.length;n++){
        if(cmult(vec[n],cconj(vec[n])).r>tol**2) return false
    }
    return true
}
function cloneMat(mat){
    var ans = []
    for(let n=0;n<mat.length;n++){
        ans[n]=[]
    for(let m=0;m<mat[0].length;m++){
        ans[n][m]=z(mat[n][m].r,mat[n][m].i)
    }
    }
    return ans
}
function crpow(c,i){
    var a = Math.atan2(c.i,c.r)
    var r = Math.hypot(c.r,c.i)
    return z(Math.cos(a*i)*(r**i),Math.sin(a*i)*(r**i))
}

function max2(x,y){
    if(Math.abs(x)<Math.abs(y)) return y
    return x
}
function normalise(vec){
    v = cloneMat([vec])[0]
    var imag = crpow(dotp(v,v),-0.5)

    for(let n=0;n<vec.length;n++){
        v[n] = cmult(imag,v[n])
    }
    return v
}
function cap1(x){
    if(Math.abs(x)>1){
        return Math.sign(x)
    }
    return x
}
function cconj(c){
    return z(c.r,-c.i)
}
function vmag(vec){
    var sum = z(0,0)
    for(let n=0;n<vec.length;n++){
        sum = cadd(sum, cmult(vec[n],cconj(vec[n])))
    }
    return crpow(sum,0.5)
}
function dotp(v1,v2){
    ans = z(0,0)
if(v1[0][0] ^ v2[0][0]) throw "Can't submit both a row & column vector"
    if(v1[0][0] && v2[0][0]){
        for(let n=0;n<Math.max(v1.length,v2.length);n++){
        ans = cadd(ans,cmult(v1[n][0],cconj(v2[n][0])))
        }
    }else{
    for(let n=0;n<Math.max(v1.length,v2.length);n++){
        ans = cadd(ans,cmult(v1[n],cconj(v2[n])))
    }
}
return ans
}
function addRows(r1,r2){
    var ans = []
    for(let n=0;n<r1.length;n++){
        ans[n]=cadd(r1[n],r2[n])
    }
    return ans
}
function subRows(r1,r2){
    var ans = []
    for(let n=0;n<r1.length;n++){
        ans[n]=csub(r1[n],r2[n])
    }
    return ans
}
function cmultRow(c,r){
    var ans = []
    for(let n=0;n<r.length;n++){
        ans[n]=cmult(c,r[n])
    }
    return ans
}

function z(real,imag){
    return new ComplexNumber(real,imag)
}
function cmult(x,y){
    return z(x.r*y.r - x.i*y.i , x.r*y.i + x.i*y.r)
}
function cadd(x,y){
    return z(x.r+y.r,x.i+y.i)
}
function csub(x,y){
    return z(x.r-y.r,x.i-y.i)
}
function mmult(m1,m2){
    ans = []
var n
var m
var o
    for(n=0;n<m2[0].length;n++){
            for(m=0;m<m1.length;m++){
                if(ans[m]==undefined)ans[m]=[]
                var sum = z(0,0)
        for(o=0;o<m2.length;o++){
//console.log(cmult(m1[m][o],m2[o][n]).r)
            sum = cadd(sum,cmult(m1[m][o],m2[o][n]))
        }
                ans[m][n]=sum
    }
    }
    return ans
}

function swapRows(mat,r1,r2){
for(let n = 0;n<mat[0].length;n++){
    var tmp = z(mat[r1][n].r,mat[r1][n].i)
    mat[r1][n]=z(mat[r2][n].r,mat[r2][n].i)
    mat[r2][n]=z(tmp.r,tmp.i)
}
}
function swapCols(mat,c1,c2){
for(let n = 0;n<mat.length;n++){
    var tmp = z(mat[n][c1].r,mat[n][c1].i)
    mat[n][c1]=z(mat[n][c2].r,mat[n][c2].i)
    mat[n][c2]=z(tmp.r,tmp.i)
}
}



function I(size){
var n
var m
    var ans = []
    for(n=0;n<size;n++){
ans.push([])
        for(m=0;m<size;m++){
//console.log(n,ans)
            if(n==m){ans[n].push(z(1,0))}else{
                        
            ans[n].push(z(0,0))
                }
        }
    }
    return ans
}

function cinv(c){
return cmult(z(c.r,-c.i),z(1/Math.hypot(c.r,c.i)**2,0))
}
function cround(x,m){
    return z(Math.round(x.r*(10**m))/(10**m),Math.round(x.i*(10**m))/(10**m))
}
function idmin(I,matr){
   var mat = cloneMat(matr)
    var ans = []
    for(let n=0;n<mat.length;n++){
        ans[n]=[]
        for(let m=0;m<mat.length;m++){
            if(n==m){
                ans[n][m]=csub(mat[n][m],z(I,0))
                continue
            }
            ans[n][m]=mat[n][m]
        }
    }
    return ans
}
function derPoly(poly){
    ans=[]
    for(let n=0;n<poly.length-1;n++){
        ans[n]=cmult(z((n+1),0),poly[n+1])
        
    }
    return ans
}

function evalPoly(poly,x){
    var X = z(1,0)
    var sum = z(0,0)
    for(let n=0;n<poly.length;n++){
        sum = cadd(cmult(X,poly[n]),sum)
        X=cmult(X,x)
    }
    return sum
}