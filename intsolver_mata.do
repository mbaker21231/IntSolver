*! Intsolver version 1.0.0 mjbaker
*! Intsolver version 1.0.1 mjbaker
*! Updated to include "usols" function 27Aug2014
	clear all
	mata:
	mata set matastrict on

	struct int_problem {
		pointer (real scalar function) scalar grad_I
		pointer (real scalar function) scalar jac_I
		pointer (real scalar function) scalar grad
		pointer (real scalar function) scalar jac
		real matrix Init_Ints, Ints, Pnts, Pandian_p1, Pandian_p2
		transmorphic addinfo
		real scalar args,digits,mbwidth,maxit,tol,useaddinfo,
				  tolsols,heuristic
	}

	void int_prob_f_Iform(struct int_problem PI,pointer (real scalar function) grad_I) PI.grad_I=grad_I

	void int_prob_jac_Iform(struct int_problem PI,pointer (real scalar function) jac_I) PI.jac_I=jac_I

	void int_prob_f(struct int_problem PI,pointer (real scalar function) grad) PI.grad=grad

	void int_prob_jac(struct int_problem PI,pointer (real scalar function) jac) PI.jac=jac

	void int_prob_args(struct int_problem PI,real scalar args) PI.args=args

	void int_prob_args_view(struct int_problem PI) PI.args

	void int_prob_mbwidth(struct int_problem PI, real scalar mbwidth) PI.mbwidth=mbwidth

	void int_prob_maxit(struct int_problem PI,real scalar maxit) PI.maxit=maxit

	void int_prob_digits(struct int_problem PI,real scalar digits) PI.digits=digits

	void int_prob_tol(struct int_problem PI,real scalar tol) PI.tol=tol

	void int_prob_tolsols(struct int_problem PI,real scalar tolsols) PI.tolsols=tolsols
	
	void int_prob_addinfo(struct int_problem PI,transmorphic addinfo) {
		PI.addinfo=addinfo
		PI.useaddinfo=1
	}
	
	void int_prob_method(struct int_problem PI, string scalar type) {
		if (type=="newton") PI.heuristic=1
	}
	
	real matrix int_prob_pts_vals(struct int_problem PI) return(PI.Pnts)
	real matrix int_prob_pandian1(struct int_problem PI) return(PI.Pandian_p1)
	real matrix int_prob_pandian2(struct int_problem PI) return(PI.Pandian_p2)

/******************************/

	void int_prob_ival(struct int_problem PI,real matrix Ints) {
		if (rows(Ints)==1) PI.Ints=Ints
		else PI.Ints=rowshape(Ints,1)
		PI.Init_Ints=PI.Ints
	}

	void int_prob_init_pts(struct int_problem PI, real matrix Pnts) {
		PI.Pnts=Pnts
	}
	
	void int_assign(struct int_problem PI, real matrix Ints) PI.Ints=Ints

	real matrix int_prob_ints_vals(struct int_problem PI) return(PI.Ints)

	struct int_problem int_prob_init()
	{
		struct int_problem scalar PI
		PI.digits=1e-10
		PI.mbwidth=1e-4
		PI.tol=1e-8
		PI.tolsols=1e-2
		PI.maxit=50
		PI.useaddinfo=0
		PI.heuristic=0
		return(PI)
	}

	real matrix make_fullgrad(struct int_problem PI,real matrix Ps)
	{
		real scalar i
		real matrix Rowgrad
	
		Rowgrad=J(rows(Ps),0,.)
		for (i=1;i<=cols(Ps);i++) Rowgrad=Rowgrad,int_evaluate_grad(PI,Ps,i)
		return(Rowgrad)
	}
	
	real matrix make_fulljac(struct int_problem PI,real matrix Ps)
	{
		real scalar i,j
		real matrix Rowjac

		Rowjac=J(rows(Ps),0,.)
		for (i=1;i<=cols(Ps);i++) {
			for (j=1;j<=cols(Ps);j++) {
				Rowjac=Rowjac,int_evaluate_jac(PI,Ps,i,j)
			}
		}
		return(Rowjac)
	}

	real matrix make_fullgradI(struct int_problem PI,real matrix Ps)
	{
		real scalar i
		real matrix Rowgrad
	
		Rowgrad=J(rows(Ps),0,.)
		for (i=1;i<=cols(Ps)/2;i++) Rowgrad=Rowgrad,int_evaluate_Igrad(PI,Ps,i)
		return(Rowgrad)
	}

	real matrix make_fulljacI(struct int_problem PI,real matrix Ps)
	{
		real scalar i,j
		real matrix Rowjac
	
		Rowjac=J(rows(Ps),0,.)
		for (i=1;i<=cols(Ps)/2;i++) {
			for (j=1;j<=cols(Ps)/2;j++) {
				Rowjac=Rowjac,int_evaluate_Ijac(PI,Ps,i,j)
			}
		}
		return(Rowjac)
	}
	
	real matrix int_evaluate_Igrad(struct int_problem PI,real matrix P,real scalar i)
	{
		if (PI.useaddinfo==0) return((*PI.grad_I)(P,i,PI.digits))
		else return((*PI.grad_I)(P,i,PI.digits, PI.addinfo))
	}

	real matrix int_evaluate_Ijac(struct int_problem PI,real matrix P,real scalar i,real scalar j)
	{
		if (PI.useaddinfo==0) return((*PI.jac_I)(P,i,j,PI.digits))
		else return((*PI.jac_I)(P,i,j,PI.digits, PI.addinfo))
	}

	real matrix int_evaluate_grad(struct int_problem PI,real matrix P,real scalar i)
	{
		if (PI.useaddinfo==0) return((*PI.grad)(P,i))
		else return((*PI.grad)(P,i,PI.addinfo))
	}
	real matrix int_evaluate_jac(struct int_problem PI,real matrix P,real scalar i,real scalar j)
	{
		if (PI.useaddinfo==0) return((*PI.jac)(P,i,j))
		else return((*PI.jac)(P,i,j,PI.addinfo))
	}

	void gsrow_intsolve(struct int_problem PI,real scalar i,real scalar j)
	{
		real matrix M,fm,Mj,fij,test,P1,P2,P3,P4,M1,M2,M3,M4,Mj1,Mj2,Mj3,Mj4,fm1,fm2,fm3,fm4,
	            fij1,fij2,fij3,fij4,dfij,dfij1,dfij2,dfij3,dfij4,Pn,Xc,fik,rhs,Pnew,P
		real scalar k,digits,n

		digits=PI.digits		
		n=PI.args
		P=PI.Ints

		P=uniqrows(P)
		M=int_mid(P,digits)
		fm=int_evaluate_Igrad(PI,M,i)
		Mj=M[.,2*j-1::2*j]
		fij=int_evaluate_Ijac(PI,P,i,j)
		if (all(fij:==0)) {
			return	/* Get out of here if there is no potential for progress */
		}

		test=sign(fij[.,1]:*fij[.,2])
		P1=select(P,test:==1)
		P2=select(P,fij[,1]:==0)
		P3=select(P,fij[,2]:==0)
	
	/* How many of these things are both? */

		P4=select(P,test:==-1)
		M1=select(M,test:==1)
		M2=select(M,fij[,1]:==0)
		M3=select(M,fij[,2]:==0)
		M4=select(M,test:==-1)
		Mj1=select(Mj,test:==1)
		Mj2=select(Mj,fij[,1]:==0)
		Mj3=select(Mj,fij[,2]:==0)
		Mj4=select(Mj,test:==-1)
		fm1=select(fm,test:==1)
		fm2=select(fm,fij[,1]:==0)
		fm3=select(fm,fij[,2]:==0)
		fm4=select(fm,test:==-1)
		fij1=select(fij,test:==1)
		fij2=select(fij,fij[,1]:==0)
		fij3=select(fij,fij[,2]:==0)
		fij4=select(fij,test:==-1)
		dfij1=(r_down(1:/fij1[,2],digits),r_up(1:/fij1[,1],digits))
		dfij2=(r_down(1:/fij2[,2],digits),J(rows(fij2),1,1e+200))
		dfij3=(J(rows(fij3),1,-1e+200),r_up(1:/fij3[,1],digits))
		dfij4=(J(rows(fij4),1,-1e+200),r_up(1:/fij4[,1],digits)),(r_down(1:/fij4[,2],digits),J(rows(fij4),1,1e+200))
		P4=colshape(J(1,2,1)#P4,2*n)
		M4=colshape(J(1,2,1)#M4,2*n)
		Mj4=colshape(J(1,2,1)#Mj4,2)
		fm4=colshape(J(1,2,1)#fm4,2)
		fij4=colshape(J(1,2,1)#fij4,2)
		dfij4=colshape(dfij4,2)
		Pn=P1\P2\P3\P4
		M=M1\M2\M3\M4
		Mj=Mj1\Mj2\Mj3\Mj4
		fm=fm1\fm2\fm3\fm4
		fij=fij1\fij2\fij3\fij4
		dfij=dfij1\dfij2\dfij3\dfij4
		Xc=int_sub(Pn,M,digits)
		fik=J(rows(Pn),2,0)	

		for (k=1;k<=n;k++) {
			if (k!=j) fik=int_add(int_mult(int_evaluate_Ijac(PI,Pn,i,k),Xc[.,2*k-1::2*k],digits),fik,digits)
		}

		rhs=int_sub(Mj,int_mult(dfij,int_add(fm,fik,digits),digits),digits)
		Pnew=Pn
		Pnew[.,2*j-1::2*j]=int_int(Pn[.,2*j-1::2*j],rhs)
		PI.Ints=uniqrows(Pnew)
	}
	
	void gsrowcycler_I(struct int_problem PI)
	{
		real matrix Ptest
		real scalar i,j,n
	
		n=PI.args
		do {
			for (i=1;i<=n;i++) {
				for (j=1;j<=n;j++) {
					PI.Ints=select(PI.Ints,rowmissing(PI.Ints):==0)
					if (rows(PI.Ints)!=0) {
						gsrow_intsolve(PI,i,j)
					}
				}
			}

			Ptest=PI.Ints						
			PI.Ints=select(PI.Ints,rowmissing(PI.Ints):==0)
	   } while (PI.Ints!=Ptest & rows(PI.Ints)!=0)
	}
	
	void int_solve(struct int_problem PI)
	{
		real matrix P,Test,rc,dc,Pd
		real scalar tolcrit 
	
		Pd=J(0,2*PI.args,.)
		rc=1..PI.args
		do {
			gsrowcycler_I(PI)
			PI.Ints=select(PI.Ints,rowmissing(PI.Ints):==0)
			if (rows(PI.Ints)!=0) {
				PI.Ints=uniqrows(PI.Ints)
				PI.Ints=int_mince(PI.Ints,rc[1],2,PI.digits)
				Test=int_evaluate_Igrad(PI,PI.Ints,rc[1])
				PI.Ints=select(PI.Ints,(Test[,1]:<=0):*(Test[,2]:>=0))
				tolcrit=rowmax(int_widths(PI.Ints))
				dc=tolcrit:<PI.mbwidth
				Pd=Pd\select(PI.Ints,dc)
				PI.Ints=select(PI.Ints,1:-dc)
				rc=rc[cols(rc)],rc[1::cols(rc)-1]
				PI.Ints=uniqrows(PI.Ints)
			}

		} while (rows(PI.Ints)>0 & max(tolcrit)>PI.mbwidth)
		PI.Ints=uniqrows(Pd\PI.Ints)
	}

	/* Functions that we need to do the rest of the stuff */
	/* First, program gradient and other functions */

	void int_newton_iter(struct int_problem PI)
	{
		real scalar it,n,i,cha1,cha2,cha
		real matrix Bds,po,pn,lbs,ubs,AA,Am,f,Amf,final

		if (PI.heuristic==0 & rows(PI.Ints)==0)  return
		else if (PI.heuristic==0) PI.Pnts=mid(PI.Ints)

		Bds=PI.Init_Ints

		po=PI.Pnts
		it=0
		n=PI.args
		lbs=J(1,0,.)
		ubs=J(1,0,.)
		for (i=1;i<=2*n;i=i+2) {
			lbs=lbs,Bds[,i]
			ubs=ubs,Bds[,i+1]
		}

		do {
			AA=make_fulljac(PI,po)
			Am=rm_newtinv(AA,PI.maxit,PI.tol)
			f=make_fullgrad(PI,po)
			Amf=rm_matvecmult(Am,f)
			pn=po:-Amf
			cha1=rowsum((pn:-po):^2)
			cha2=rowsum(po:^2)
			cha=max(cha1:/cha2)
			if (hasmissing(cha)) cha=max(cha1)
			po=pn
			it=it+1
			po=select(po,rowsum(po:>=lbs):==n)
			po=select(po,rowsum(po:<=ubs):==n)
		} while (cha>PI.tol & it<=PI.maxit & rows(po)>0)
		
		if (rows(po)>1) po=usols(po,PI.tolsols)
		PI.Pnts=po
	}

	void int_pandian(struct int_problem PI)
	{
		real scalar i
		real matrix X,eta,IJac,L,M,Ihat,rho,ind,
			Selmat,lhs,A,Am,P,C,P_new
	
		Ihat=PI.Ints:+J(rows(PI.Ints),cols(PI.Ints)/2,(-PI.tol,PI.tol))

		X=mid(PI.Ints)
		rho=J(rows(X),0,.)
		for (i=1;i<cols(Ihat);i=i+2) rho=rho,radius(Ihat[,i::i+1])
		IJac=make_fulljacI(PI,Ihat)
		L=J(rows(IJac),0,.)
		M=J(rows(IJac),0,.)
		for (i=1;i<cols(IJac);i=i+2) {
			L=L,IJac[,i]
			M=M,IJac[,i+1]
		}
	
		A=(M+L)/2
		Am=rm_newtinv(A,PI.maxit,PI.tol)
		P=abs(Am)
		eta=abs(rm_matvecmult(Am,make_fullgrad(PI,X)))
		L=abs(A-L)
		M=abs(A-M)
		C=J(rows(L),0,.)
		for (i=1;i<=cols(L);i++) C=C,rowmax((L[,i],M[,i]))
		P=rm_matmult(P,C)
		lhs=eta:+rm_matvecmult(P,rho)
		ind=rowsum(lhs:<=rho):==cols(X)
		PI.Pandian_p1=ind
	
	/* Second part of test */

		Selmat=rowshape(lowertriangle(J(sqrt(cols(L)),sqrt(cols(L)),1)),1)
		Selmat=J(rows(P),1,Selmat)
		P_new=P:*Selmat
		lhs=rm_matvecmult(P_new,rho)
		ind=rowsum(lhs:<rho):==cols(X)
		PI.Pandian_p2=ind
	}

	real matrix usols(real matrix p,real scalar coltol)
	{
		real matrix pdone,test
		real scalar i

		if (rows(p)<=1) return(p)
		pdone=p[1,.]
		for (i=2;i<=rows(p);i++) {
			test=rowsum(abs(pdone:-J(rows(pdone),1,p[i,.])))/cols(pdone)
			if (all(test:>coltol)) pdone=pdone \ p[i,.]
		}
		return(pdone)
	}
	
	mata mlib create lintsolver, dir(PERSONAL) replace

	mata mlib add lintsolver *()

	mata mlib index

end
