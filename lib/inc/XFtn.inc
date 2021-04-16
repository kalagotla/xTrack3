c
c	Function declarations.
c
	integer*8 XFtnInit, XFtnGC, XFtnEvent,
     &            XFtnImgCreate, XFtnMakeDrawable, XFtnWRead,
     &            XFtnColorWindow, XFtnGLWindow, XFtnReport
	integer XFtnNumEvents, XFtnInitDial, XFtnGetDials,
     &          XFtnGetEvent
c
c	Event names.  Use in `type' field in XEvent structure.
c
	integer KeyPress, KeyRelease, ButtonPress, ButtonRelease,
     &	  MotionNotify, EnterNotify, LeaveNotify, FocusIn, FocusOut,
     &	  KeymapNotify, Expose, GraphicsExpose, NoExpose, CreateNotify,
     &	  DestroyNotify, UnmapNotify, MapNotify, MapRequest, 
     &	  ReparentNotify, ConfigureNotify, ConfigureRequest, 
     &	  GravityNotify, ResizeRequest, CirculateNotify, 
     &	  CirculateRequest, PropertyNotify, SelectionClear,
     &	  SelectionRequest, SelectionNotify, ColormapNotify,
     &	  ClientMessage, MappingNotify, LASTEvent

	parameter ( KeyPress		=  2 )
	parameter ( KeyRelease		=  3 )
	parameter ( ButtonPress		=  4 )
	parameter ( ButtonRelease	=  5 )
	parameter ( MotionNotify	=  6 )
	parameter ( EnterNotify		=  7 )
	parameter ( LeaveNotify		=  8 )
	parameter ( FocusIn		=  9 )
	parameter ( FocusOut		= 10 )
	parameter ( KeymapNotify	= 11 )
	parameter ( Expose		= 12 )
	parameter ( GraphicsExpose	= 13 )
	parameter ( NoExpose	 	= 14 )
	parameter ( CreateNotify	= 16 )
	parameter ( DestroyNotify	= 17 )
	parameter ( UnmapNotify		= 18 )
	parameter ( MapNotify		= 19 )
	parameter ( MapRequest		= 20 )
	parameter ( ReparentNotify	= 21 )
	parameter ( ConfigureNotify	= 22 )
	parameter ( ConfigureRequest	= 23 )
	parameter ( GravityNotify	= 24 )
	parameter ( ResizeRequest	= 25 )
	parameter ( CirculateNotify	= 26 )
	parameter ( CirculateRequest	= 27 )
	parameter ( PropertyNotify	= 28 )
	parameter ( SelectionClear	= 29 )
	parameter ( SelectionRequest	= 30 )
	parameter ( SelectionNotify	= 31 )
	parameter ( ColormapNotify	= 32 )
	parameter ( ClientMessage	= 33 )
	parameter ( MappingNotify	= 34 )
	parameter ( LASTEvent		= 35 )
c
c	Button masks.  Used as modifiers to Grab(Button,Key),
c	results of QueryPointer.  Not to be confused with button names below.
c
	integer Button1Mask, Button2Mask, Button3Mask, Button4Mask,
     &	    Button5Mask

	parameter ( Button1Mask =  256 )
	parameter ( Button2Mask =  512 )
	parameter ( Button3Mask = 1024 )
	parameter ( Button4Mask = 2048 )
	parameter ( Button5Mask = 4096 )
c
c	Button names.  Used as arguments in GrabButton and as detail
c	in Button(Press,Release) events.  Not to be confused with 
c	button masks above.
c
	integer Button1, Button2, Button3, Button4, Button5

	parameter ( Button1	= 1 )
	parameter ( Button2	= 2 )
	parameter ( Button3	= 3 )
	parameter ( Button4	= 4 )
	parameter ( Button5	= 5 )
c
c	Graphics functions
c
	integer GXclear, GXand, GXandReverse, GXcopy, GXandInverted,
     &	    GXnoop, GXxor, GXor, GXnor, GXequiv, GXinvert,
     &	    GXorReverse, GXcopyInverted, GXorInverted, GXnand, GXset

	parameter ( GXclear		=  0 )
c						0
	parameter ( GXand		=  1 )
c						src AND dst
	parameter ( GXandReverse	=  2 )
c 						src AND NOT dst
	parameter ( GXcopy		=  3 )
c 						src
	parameter ( GXandInverted	=  4 )
c 						NOT src AND dst
	parameter ( GXnoop		=  5 )
c 						dst
	parameter ( GXxor		=  6 )
c 						src XOR dst
	parameter ( GXor		=  7 )
c 						src OR dst
	parameter ( GXnor		=  8 )
c 						NOT src AND NOT dst
	parameter ( GXequiv		=  9 )
c 						NOT src XOR dst
	parameter ( GXinvert		= 10 )
c 						NOT dst
	parameter ( GXorReverse		= 11 )
c 						src OR NOT dst
	parameter ( GXcopyInverted	= 12 )
c 						NOT src
	parameter ( GXorInverted	= 13 )
c 						NOT src OR dst
	parameter ( GXnand		= 14 )
c 						NOT src OR NOT dst
	parameter ( GXset		= 15 )
c 						1
c
c	Matrix types.
c
	integer XFDIObjectMatrix, XFDINormalMatrix, XFDILightMatrix

	parameter ( XFDIObjectMatrix = 0 )
	parameter ( XFDINormalMatrix = 1 )
	parameter ( XFDILightMatrix  = 2 )
c
c	Matrix flags.
c
	integer XFDIMatrix4X4, XFDIMatrix4X3, XFDIMatrix3X3,
     &          XFDIMatrix2D, XFDIMatrixXRotate, XFDIMatrixYRotate, 
     &          XFDIMatrixZRotate, XFDIMatrixScale, XFDIMatrixTranslate,
     &          XFDIMatrixIdentity

	parameter ( XFDIMatrix4X4       = 0 )
	parameter ( XFDIMatrix4X3       = 1 )
	parameter ( XFDIMatrix3X3       = 2 )
	parameter ( XFDIMatrix2D        = 3 )
	parameter ( XFDIMatrixXRotate   = 4 )
	parameter ( XFDIMatrixYRotate   = 5 )
	parameter ( XFDIMatrixZRotate   = 6 )
	parameter ( XFDIMatrixScale     = 7 )
	parameter ( XFDIMatrixTranslate = 8 )
	parameter ( XFDIMatrixIdentity  = 9 )
c
	integer XFDIOff, XFDIOn, XFDINotThere

	parameter ( XFDIOff      = 0 )
	parameter ( XFDIOn       = 1 )
	parameter ( XFDINotThere = 2 )
