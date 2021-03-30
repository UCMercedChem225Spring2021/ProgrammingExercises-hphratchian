      Program PIB
!
!     To compile this program using pgfortran, run the command:
!           pgfortran -llapack -lblas -o pib.exe pib.f03
!
!     This program carries out a variational calculation for a quantum
!     particle-in-box modified by a linear potential. Specifically, the
!     potential is given by:
!
!           V(x) = b x      for  0 < x < l
!           V(x) = inf      for  x <= 0  and  x >= l
!
!     where b and l are user-provided parameters.
!
!     Planck's constant and the particle mass are both taken to be 1.0.
!
!
!     The variational problem is solved in the particle-in-a-box
!     eigenfunction basis. The user provides the number of basis functions
!     to be used, NBasis; the basis set is take to be the first NBasis
!     particle-in-a-box eigenfunctions.
!
!
!     H.P. Hratchian, 2016.
!
!
!     Variable Declarations
      Implicit None
      Integer::i,NCmdLineArgs,NBasis
      Real::l,b,mass,start_time_total,end_time_total,start_time_local,  &
        end_time_local
      Real,Dimension(:),Allocatable::HEVals
      Real,Dimension(:,:),Allocatable::TMat,VMat,HMat,HEVecs
      Logical::Read_CmdLine_Value,Read_CmdLine_l,Read_CmdLine_b,  &
        Read_CmdLine_Mass,Read_CmdLine_NBasis,Print_Arrays
      Character(Len=1024)::cmd_buffer
!
!     Format Statements
 1000 Format(1x,'What is the value of l (box length)?')
 1010 Format(1x,'What is the value of b (the slope of the linear ',  &
        'potential)?')
 1020 Format(1x,'How many basis functions are included?')
 2000 Format(/,1x,'The reduced Planck''s constant is set to 1.')
 2010 Format(/,1x,'Box Length      = ',F10.3,/,  &
        1x,'Potential Slope = ',F10.3,/,  &
        1x,'Mass            = ',F10.3,/,  &
        1x,'NBasis          = ',I12,/)
 8000 Format(1x,A,': ',F10.1,' s')
 9000 Format(/,1x,'Confused command line reading...',  &
        'found switch but unsure what to read next.')
 9010 Format(/,1x,'Found unknown command line argument: ',A)
!
!
!     Begin by setting defaults for input parameters and then processing
!     the command line arguments to overide these parameters as directed by
!     the user.
!
      Call CPU_Time(start_time_total)
      l = Float(1)
      b = Float(0)
      mass = Float(1)
      NBasis = 5
      Print_Arrays = .True.
      NCmdLineArgs = command_argument_count()
      Read_CmdLine_Value = .False.
      Read_CmdLine_l     = .False.
      Read_CmdLine_b     = .False.
      Do i = 1,NCmdLineArgs
        Call Get_Command_Argument(i,cmd_buffer)
        cmd_buffer = AdjustL(cmd_buffer)
        If(Read_CmdLine_Value) then
          Read_CmdLine_Value = .False.
          If(Read_CmdLine_l) then
            Read(cmd_buffer,*) l
            Read_CmdLine_l = .False.
          elseif(Read_CmdLine_b) then
            Read(cmd_buffer,*) b
            Read_CmdLine_b = .False.
          elseif(Read_CmdLine_mass) then
            Read(cmd_buffer,*) mass
            Read_CmdLine_mass = .False.
          elseif(Read_CmdLine_NBasis) then
            Read(cmd_buffer,*) NBasis
            Read_CmdLine_NBasis = .False.
          else
            Write(*,9000)
            STOP
          endIf
        else
          Select Case(cmd_buffer)
            Case("-noprintarrays")
              Print_Arrays = .False.
            Case("-printarrays")
              Print_Arrays = .True.
            Case("-l","-length")
              Read_CmdLine_Value = .True.
              Read_CmdLine_l = .True.
            Case("-b","-slope")
              Read_CmdLine_Value = .True.
              Read_CmdLine_b = .True.
            Case("-m","-mass")
              Read_CmdLine_Value = .True.
              Read_CmdLine_mass = .True.
            Case("-nbasis")
              Read_CmdLine_Value = .True.
              Read_CmdLine_NBasis = .True.
            Case Default
              Write(*,9010) Trim(cmd_buffer)
              STOP
          End Select
        endIf
      EndDo
      Write(*,2000)
      Write(*,2010) l,b,mass,NBasis
!
!     Allocate memory for the kinetic, potential, Hamiltonian, and
!     Hamiltonian eigen-values/vectors arrays. Then, evaluate the integrals
!     and fill the matrices.
!
      Allocate(TMat(NBasis,NBasis),VMat(NBasis,NBasis),  &
        HMat(NBasis,NBasis))
      Allocate(HEVals(NBasis),HEVecs(NBasis,NBasis))
      Call CPU_Time(start_time_local)
      Call Fill_PIB_TMat(NBasis,l,mass,TMat)
      Call CPU_Time(end_time_local)
      Write(*,8000)'Time for TMat formation',end_time_local - start_time_local
      Call CPU_Time(start_time_local)
      Call Fill_PIB_VMat(NBasis,l,b,VMat)
      Call CPU_Time(end_time_local)
      Write(*,8000)'Time for VMat formation',end_time_local - start_time_local
      Call CPU_Time(start_time_local)
      HMat = TMat + VMat
      Call CPU_Time(end_time_local)
      Write(*,8000)'Time for (TMat + VMat) formation',end_time_local - start_time_local
      If(Print_Arrays) then
        Call Print_Matrix_Full_Real(6,TMat,'T:',NBasis,NBasis)
        Call Print_Matrix_Full_Real(6,VMat,'V:',NBasis,NBasis)
        Call Print_Matrix_Full_Real(6,HMat,'H:',NBasis,NBasis)
      endIf
!
!     Diagonalize the Hamiltonian. The eigen-values/vectors will be loaded
!     into arrays HEVals and HEVecs.
!
      Call CPU_Time(start_time_local)
      Call Matrix_Diagonalize(NBasis,HMat,HEVecs,HEVals)
      Call CPU_Time(end_time_local)
      Write(*,8000)'Time for H diagonalization',end_time_local - start_time_local
      If(Print_Arrays) then
        Call Print_Matrix_Full_Real(6,HEVals,  &
          'Hamiltonian Eigen-Values:',NBasis,1)
        Call Print_Matrix_Full_Real(6,HEVecs,  &
          'Hamiltonian Eigen-Vectors:',NBasis,NBasis)
      endIf
!
!     Stop the total-job clock and report job time.
!
      Call CPU_Time(end_time_total)
      Write(*,8000) 'Total Job Time',end_time_total - start_time_total
      End Program PIB


!PROCEDURE Fill_PIB_TMat
      Subroutine Fill_PIB_TMat(NBasis,l,mass,TMat)
!
!     This subroutine fills the kinetic energy matrix for particle-in-a-box
!     basis functions running from n=1 through n=NBasis (where n is the
!     quantum number). It is assumed that the "box" has a length l and that
!     the internal units are based on h=1.
!
!
!     H.P. Hratchian, 2016.
!
!
!     Variable Declarations
      Implicit None
      Integer,Intent(In)::NBasis
      Real,Intent(In)::l,mass
      Real,Dimension(NBasis,NBasis),Intent(InOut)::TMat
!
      Integer::i
      Real::pi,Prefactor
!
!
!     Initialize TMat to 0.0 and fill the diagonal.
!
      pi=float(4)*ATan(float(1))
      Prefactor = (pi**2)/(float(2)*mass*l**2)
      TMat = Float(0)
      Do i = 1,NBasis
        TMat(i,i) = Prefactor*Float(i)**2
      EndDo
!
      End Subroutine Fill_PIB_TMat


!PROCEDURE Fill_PIB_VMat
      Subroutine Fill_PIB_VMat(NBasis,l,b,VMat)
!
!     This subroutine fills the potential energy matrix for
!     particle-in-a-box basis functions running from n=1 through n=NBasis
!     (where n is the quantum number). It is assumed that the "box" has a
!     length l and that the internal units are based on h=1.
!
!
!     H.P. Hratchian, 2016.
!
!
!     Variable Declarations
      Implicit None
      Integer,Intent(In)::NBasis
      Real,Intent(In)::l,b
      Real,Dimension(NBasis,NBasis),Intent(InOut)::VMat
!
      Integer::i,j,IntTemp
      Real::One,Pi,Prefactor,RealTemp
!
!
!     Initialize VMat to 0.0 and fill the diagonal.
!
      One = Float(1)
      Pi = Float(4)*ATan(One)
      Prefactor = b*l/Pi**2
      Do i = 1,NBasis
        Do j = 1,NBasis
          If(i.ne.j) then
            IntTemp = (i-j)
            RealTemp = Float(IntTemp)*Pi
            RealTemp = Cos(RealTemp)-One
            IntTemp = IntTemp**2
            VMat(i,j) = One/Float(IntTemp)
            IntTemp = (i+j)**2
            VMat(i,j) = VMat(i,j)-One/Float(IntTemp)
            VMat(i,j) = Prefactor*RealTemp*VMat(i,j)
          else
            VMat(i,j) = b*l/Float(2)
          endIf
        EndDo
      EndDo
!
      End Subroutine Fill_PIB_VMat


!
!PROCEDURE Print_Matrix_Full_Real
      Subroutine Print_Matrix_Full_Real(IOut,A,Header,M,N)
!
!     This subroutine prints a real matrix that is fully dimension
!     - i.e., not stored in packed form.  The unit number for the output
!     file is IOut.  A is the matrix, which is dimensioned (m,n).  Header
!     is a character string that is printed at the top of the matrix dump.
!
!
!     Variable Declarations
!
      implicit none
      integer,intent(in)::IOut,M,N
      real,dimension(M,N),intent(in)::A
      character(len=*),intent(in)::Header
!
!     Local variables
      integer,parameter::NColumns=5
      integer::i,j,IFirst,ILast
!
 1000 Format(1x,A)
 2000 Format(5x,5(7x,I7))
 2010 Format(1x,I7,5F14.6)
!
      write(IOut,1000) TRIM(Header)
      do IFirst = 1,N,NColumns
        ILast = Min(IFirst+NColumns-1,N)
        write(IOut,2000) (i,i=IFirst,ILast)
        do i = 1,M
          write(IOut,2010) i,(A(i,j),j=IFirst,ILast)
        enddo
      enddo
!
      Return
      End Subroutine Print_Matrix_Full_Real


!
!     PROCEDURE Matrix_Diagonalize
      Subroutine Matrix_Diagonalize(N,A,A_EVecs,A_EVals)
!
!     This routine carries out matrix diagonalization of a symmetric square
!     matrix with leading dimension N. It is assumed that A is sent here as
!     an (N x N) array.
!
      Implicit None
      Integer,Intent(In)::N
      Real,Dimension(N,N),Intent(In)::A
      Real,Dimension(N,N),Intent(Out)::A_EVecs
      Real,Dimension(N),Intent(Out)::A_EVals
!
      Integer::i,j,k,IError=1
      Real::time1,time2
      Real,Dimension(:),Allocatable::A_Symm,Temp_Vector
!
!     Do the work...
!
      Call CPU_Time(time1)
      Allocate(A_Symm((N*(N+1))/2),Temp_Vector(3*N))
      k = 0
      Do i = 1,Size(A,1)
        Do j = 1,i
          k = k+1
          A_Symm(k) = A(i,j)
        EndDo
      EndDo
      Call CPU_Time(time2)
      Write(*,*)' time2-time1: ',time2-time1
      Call CPU_Time(time1)
      Call SSPEV('V','U',N,A_Symm,A_EVals,A_EVecs,N, &
        Temp_Vector,IError)
      Call CPU_Time(time2)
      Write(*,*)' DSPEV Time: ',time2-time1
      If(IError.ne.0) Write(*,'(1X,A,I4)')  &
        'DIAGONALIZATION FAILED: IError =',IError
      DeAllocate(A_Symm)
!
      Return
      End Subroutine Matrix_Diagonalize
