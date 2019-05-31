#CDEBUGFLAGS = -g
CDEBUGFLAGS = -O2
CFLAGS = ${CDEBUGFLAGS}
CC = cc

#----------------------------------------------------------*

NM_HDRS_ALL = nr.h 
NM_CSRS_ALL = nr.c
NM_OBJS_ALL = nr.o

SP_HDRS = sa.h sp.h sn.h sgy.h swb.h ${NM_HDRS_ALL}
SP_CSRS =      sp.c sn.c sgy.c swb.c ${NM_CSRS_ALL}
SP_OBJS =      sp.o sn.o sgy.o swb.o ${NM_OBJS_ALL}

sp: ${SP_OBJS}
	${CC} ${CFLAGS} ${SP_OBJS} \
-lm -o sp.ex

sp.o: sp.c sp.h sa.h
	${CC} ${CFLAGS} -c sp.c
sn.o: sn.c sn.h sa.h sgy.h
	${CC} ${CFLAGS} -c sn.c
sgy.o: sgy.c sn.h sa.h sgy.h
	${CC} ${CFLAGS} -c sgy.c
swb.o: sn.c sn.h sa.h swb.h
	${CC} ${CFLAGS} -c swb.c
slr.o: sn.c sn.h sa.h slr.h
	${CC} ${CFLAGS} -c slr.c

clean:
	/bin/rm -f *.o

.KEEP_STATE:
 
