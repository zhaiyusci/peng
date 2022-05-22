OBJDIR ?= build

.PHONY: all

all: 
	@echo "Build prerequest ... "
	@eval "cd external; ./prerequest.sh; cd ..;"
	@echo "done."
	@mkdir -p $(OBJDIR)
	@cp src/Makefile.build $(OBJDIR)/Makefile
	+@make -C $(OBJDIR)

clean:
	rm -rf $(OBJDIR)
