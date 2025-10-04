package CrazyE.graphic;

import arc.*;
import arc.scene.*;
import arc.util.*;
import mindustry.ui.*;
import mindustry.ui.fragments.MinimapFragment;
import CrazyE.graphic.CEMinimapField.*;

import static mindustry.Vars.*;
import static CrazyE.graphic.CEMinimapRenderer.*;

public class CEMinimapFragment extends MinimapFragment{
    @Override
    public void build(Group parent) {
        super.build(parent);

        CEMinimapRenderer render = (CEMinimapRenderer)renderer.minimap;

        parent.fill(t -> {
            t.setFillParent(true);
            t.visible(this::shown);
            t.update(() -> t.setBounds(0, 0, Core.graphics.getWidth(), Core.graphics.getHeight()));

            int rowCounter = 0;

            t.marginRight(10f).marginBottom(10f);

            t.add("Field type: ").style(Styles.outlineLabel);
            t.row();
            for(FieldType ft : FieldType.values()){
                if(rowCounter >= 3){
                    t.row();
                    rowCounter %= 3;
                }

                t.button(b -> {
                    b.label(ft::name);
                    b.setChecked(Structs.eq(ft, render.field));
                }, () -> render.field = ft);
                rowCounter++;
            }

            t.row();

            t.add("Map type: ").style(Styles.outlineLabel);
            t.row();
            rowCounter = 0;
            for(String s : fields.keys()){
                if(rowCounter == 3) t.row();
                t.button(b -> {
                    b.label(() -> s);
                    b.setChecked(Structs.eq(s, render.currentMap));
                }, () -> render.changed(s));
                rowCounter++;
            }
        });
    }
}
