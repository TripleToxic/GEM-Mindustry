package CrazyE.world.DynamicEnvironment;

import arc.Events;
import arc.func.Cons;
import arc.math.Mathf;
import arc.math.geom.*;
import arc.struct.Seq;
import arc.util.Log;
import arc.util.Time;
import arc.util.Tmp;
import mindustry.game.Team;
import mindustry.game.Teams;
import mindustry.gen.Unit;
import mindustry.io.SaveFileReader.*;
import mindustry.world.Tiles;

import java.io.DataInput;
import java.io.DataOutput;
import java.io.IOException;

import static java.lang.Math.*;
import static mindustry.Vars.*;
import static mindustry.game.EventType.*;
import static mindustry.io.SaveVersion.addCustomChunk;

public class Gravity implements CustomChunk{
    private static final double cons = cbrt(2),
            w1 = -(float)(cons/(2d - cons)),
            w2 = (float)(1d/(2d - cons));
    private static final float c = 25f * tilesize / 60f, maxSpeed = 0.9f*c, csquared = c * c, G = 0.1f, k = 2f, timescale = 1/128f,
    c1 = (float)(0.5d*w2),
    c2 = (float)(0.5d*(w1+w2)),
    c3 = c2,
    c4 = c1,

    d1 = (float)w2,
    d2 = (float)w1,
    d3 = (float)w2
    ;

    private int sizeX, sizeY;
    private GEMgrid PotentialField, vPotentialField, BufferField, CurrentField;
    private boolean init = false;

    String s = "";

    public Gravity(){
        Events.on(WorldLoadEvent.class, e -> {
            this.sizeX = world.width();
            this.sizeY = world.height();

            PotentialField = new GEMgrid(sizeY, sizeX);
            vPotentialField = new GEMgrid(sizeY, sizeX);
            BufferField = new GEMgrid(sizeY, sizeX);
            CurrentField = new GEMgrid(sizeY, sizeX);

            if(!init)Events.run(Trigger.beforeGameUpdate, this::update);
            if(!init)Events.on(SaveWriteEvent.class, e2 -> Log.info(s));
            init = true;
        });

        addCustomChunk("Gravity", this);
    }

    public Vec2 Accel(int y, int x){
        float Ex, Ey, B, xcomp = Tmp.v1.x, ycomp = Tmp.v1.y, xcomp2, ycomp2, dt = Time.delta*timescale;

        Ex = -0.5f*(PotentialField.get(y, x + 1).x - PotentialField.get(y, x - 1).x);
        Ey = -0.5f*(PotentialField.get(y + 1, x).x - PotentialField.get(y - 1, x).x);

        Ex -= vPotentialField.get(y, x).y;
        Ey -= vPotentialField.get(y, x).z;

        B = 0.5f*(PotentialField.get(y, x + 1).z - PotentialField.get(y, x - 1).z)
                - 0.5f*(PotentialField.get(y + 1, x).y - PotentialField.get(y - 1, x).y);

        xcomp += (Ex + k*Tmp.v1.y*B)*dt;
        ycomp += (Ey - k*Tmp.v1.x*B)*dt;

        for(int i=0; i<7; i++){
            xcomp2 = Tmp.v1.x + (Ex + 0.5f*k*(Tmp.v1.y + ycomp)*B)*dt;
            ycomp2 = Tmp.v1.y + (Ey - 0.5f*k*(Tmp.v1.x + xcomp)*B)*dt;

            xcomp = xcomp2; ycomp = ycomp2;
        }

        return Tmp.v2.set(xcomp, ycomp).scl(dt);
    }

    public void update(){
        Teams.TeamData D;
        Seq<Unit> U;
        Unit u;
        float x, y, m, w, h;
        int ix, iy;
        CurrentField.clear();
        for (Team t : Team.all) {
            D = t.data();
            if(D == null) continue;
            U = D.units;
            for(int i = 0; i < U.size; i++){
                u = U.get(i);
                x = u.x/tilesize;
                y = u.y/tilesize;
                if (x < 0.5f || x >= sizeX-0.5f || y < 0.5f || y >= sizeY-0.5f || Float.isNaN(x) || Float.isNaN(y)) continue;
                if(u.vel.len() >= maxSpeed) u.vel.setLength(maxSpeed);
                Tmp.v1.set(u.vel);

                m = G * u.mass() / Mathf.sqrt(1f - Tmp.v1.len2()/csquared);
                ix = round(x); iy = round(y);
                w = x - ix; h = y - iy;

                CurrentField.get(iy-1, ix-1).add(csquared * m, Tmp.v1.x * m * (0.5f-w)*(0.5f-h), Tmp.v1.y * m * (0.5f-w)*(0.5f-h));
                CurrentField.get(iy-1, ix).add(csquared * m, Tmp.v1.x * m * (0.5f+w)*(0.5f-h), Tmp.v1.y * m * (0.5f+w)*(0.5f-h));
                CurrentField.get(iy, ix-1).add(csquared * m, Tmp.v1.x * m * (0.5f-w)*(0.5f+h), Tmp.v1.y * m * (0.5f-w)*(0.5f+h));
                CurrentField.get(iy, ix).add(csquared * m, Tmp.v1.x * m*(0.5f+w)*(0.5f+h), Tmp.v1.y * m * (0.5f+w)*(0.5f+h));
            }
        }

        float dt = Time.delta*timescale;

        //Yoshida Integrator

        PotentialField.addmul(vPotentialField, c1*dt);
        BufferField.laplacian(PotentialField, vPotentialField).sub(CurrentField);
        vPotentialField.addmul(BufferField, d1*dt);

        PotentialField.addmul(vPotentialField, c2*dt);
        BufferField.laplacian(PotentialField, vPotentialField).sub(CurrentField);
        vPotentialField.addmul(BufferField, d2*dt);

        PotentialField.addmul(vPotentialField, c3*dt);
        BufferField.laplacian(PotentialField, vPotentialField).sub(CurrentField);
        vPotentialField.addmul(BufferField, d3*dt);

        PotentialField.addmul(vPotentialField, c4*dt);

        PotentialField.normalize();

        for (Team t : Team.all){
            D = t.data();
            if (D == null) continue;
            U = D.units;
            for (int i = 0; i < U.size; i++){
                u = U.get(i);
                x = u.x/tilesize;
                y = u.y/tilesize;
                ix = round(x); iy = round(y);
                w = x - ix; h = y - iy;

                if (x < 1.5 || x >= sizeX - 1.5 || y < 1.5 || y >= sizeY - 1.5 || Float.isNaN(x) || Float.isNaN(y)) continue;

                u.vel.add(Accel(iy - 1, ix - 1).scl((0.5f - w) * (0.5f - h)));
                u.vel.add(Accel(iy - 1, ix).scl((0.5f + w) * (0.5f - h)));
                u.vel.add(Accel(iy, ix - 1).scl((0.5f - w) * (0.5f + h)));
                u.vel.add(Accel(iy, ix).scl((0.5f + w) * (0.5f + h)));
            }
        }
    }

    @Override
    public void write(DataOutput write) throws IOException{
        if (PotentialField == null || vPotentialField == null) return;

        write.writeInt(sizeX);
        write.writeInt(sizeY);

        PotentialField.eachTile(v -> {
            try{
                write.writeFloat(v.x);
                write.writeFloat(v.y);
                write.writeFloat(v.z);
            }catch(IOException e){
                throw new RuntimeException(e);
            }
        });
        vPotentialField.eachTile(v -> {
            try{
                write.writeFloat(v.x);
                write.writeFloat(v.y);
                write.writeFloat(v.z);
            }catch(IOException e){
                throw new RuntimeException(e);
            }
        });
    }

    @Override
    public void read(DataInput read) throws IOException{
        sizeX = read.readInt();
        sizeY = read.readInt();

        PotentialField = new GEMgrid(sizeY, sizeX);
        vPotentialField = new GEMgrid(sizeY, sizeX);
        BufferField = new GEMgrid(sizeY, sizeX);
        CurrentField = new GEMgrid(sizeY, sizeX);

        PotentialField.eachTile(v -> {
            try{
                v.x = read.readFloat();
                v.y = read.readFloat();
                v.z = read.readFloat();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        });
        vPotentialField.eachTile(v -> {
            try{
                v.x = read.readFloat();
                v.y = read.readFloat();
                v.z = read.readFloat();
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        });
    }

    public static class GEMgrid {
        private final Vec3[] grid;
        private final Vec3 buffer = new Vec3();
        public int sizeY, sizeX;

        public GEMgrid(int sizeY, int sizeX) {
            grid = new Vec3[sizeX * sizeY];
            int totalSize = sizeX * sizeY;
            for (int i = 0; i < totalSize; i++){
                grid[i] = new Vec3();
            }
            this.sizeY = sizeY;
            this.sizeX = sizeX;
        }

        public final Vec3 get(int y, int x){
            return grid[x + y * sizeX];
        }

        public GEMgrid clear() {
            eachTile(v -> v.set(Vec3.Zero));
            return this;
        }

        /**
         * Note that the GEMgrid object that called this method will have its value be completely replaced.
         *
         * @param grid:  the grid that will be the values for the laplacian.
         * @param Bgrid: grid that will be used to impose the Neumann Boundary Condition: ∂(grid)/∂n = -(Bgrid)/c.
         */
        public GEMgrid laplacian(GEMgrid grid, GEMgrid Bgrid) {
            for (int i = 1; i < sizeY - 1; i++) {
                for (int j = 1; j < sizeX - 1; j++){
                    get(i, j).set(grid.get(i, j)).scl(-4f)
                            .add(grid.get(i, j - 1))
                            .add(grid.get(i, j + 1))
                            .add(grid.get(i - 1, j))
                            .add(grid.get(i + 1, j))
                            .scl(csquared);
                }
            }

            for (int i = 1; i < sizeY - 1; i++) { //X-edge
                buffer.set(Bgrid.get(i, 0)).scl(2f / c);
                get(i, 0).set(grid.get(i, 1)).scl(-4f)
                        .add(grid.get(i, 1))
                        .add(grid.get(i, 1))
                        .add(grid.get(i - 1, 0))
                        .add(grid.get(i + 1, 0))
                        .add(buffer).scl(csquared);

                buffer.set(Bgrid.get(i, sizeX - 1)).scl(2f / c);
                get(i, sizeX - 1).set(grid.get(i, sizeX - 2)).scl(-4f)
                        .add(grid.get(i, sizeX - 2))
                        .add(grid.get(i, sizeX - 2))
                        .add(grid.get(i - 1, sizeX - 1))
                        .add(grid.get(i + 1, sizeX - 1))
                        .add(buffer).scl(csquared);
            }


            for (int j = 1; j < sizeX - 1; j++) { //Y-edge
                buffer.set(Bgrid.get(0, j)).scl(2f / c);
                get(0, j).set(grid.get(1, j)).scl(-4f)
                        .add(grid.get(1, j))
                        .add(grid.get(1, j))
                        .add(grid.get(0, j - 1))
                        .add(grid.get(0, j + 1))
                        .add(buffer).scl(csquared);

                buffer.set(Bgrid.get(sizeY - 1, j)).scl(2f / c);
                get(sizeY - 1, j).set(grid.get(sizeY - 2, j)).scl(-4f)
                        .add(grid.get(sizeY - 2, j))
                        .add(grid.get(sizeY - 2, j))
                        .add(grid.get(sizeY - 1, j - 1))
                        .add(grid.get(sizeY - 1, j + 1))
                        .add(buffer).scl(csquared);
            }


            //Bottom-left
            get(0, 0).set(Bgrid.get(0, 0)).scl(1f / c).sub(grid.get(0, 0)).scl(2f)
                    .add(grid.get(0, 1)).add(grid.get(1, 0)).scl(2f * csquared);
            //Bottom-right
            get(0, sizeX - 1).set(Bgrid.get(0, sizeX - 1)).scl(1f / c).sub(grid.get(0, sizeX - 1)).scl(2f)
                    .add(grid.get(0, sizeX - 2)).add(grid.get(1, sizeX - 1)).scl(2f * csquared);
            //Top-left
            get(sizeY - 1, 0).set(Bgrid.get(sizeY - 1, 0)).scl(1f / c).sub(grid.get(sizeY - 1, 0)).scl(2f)
                    .add(grid.get(sizeY - 1, 1)).add(grid.get(sizeY - 2, 0)).scl(2f * csquared);
            //Top-right
            get(sizeY - 1, sizeX - 1).set(Bgrid.get(sizeY - 1, sizeX - 1)).scl(1f / c).sub(grid.get(sizeY - 1, sizeX - 1)).scl(2f)
                    .add(grid.get(sizeY - 1, sizeX - 2)).add(grid.get(sizeY - 2, sizeX - 1)).scl(2f * csquared);

            return this;
        }

        public GEMgrid set(GEMgrid grid) {
            for (int i = 0; i < sizeY; i++) {
                for (int j = 0; j < sizeX; j++) {
                    get(i, j).set(grid.get(i, j));
                }
            }
            return this;
        }

        public GEMgrid add(GEMgrid grid) {
            for (int i = 0; i < sizeY; i++) {
                for (int j = 0; j < sizeX; j++) {
                    get(i, j).add(grid.get(i, j));
                }
            }
            return this;
        }

        public GEMgrid sub(GEMgrid grid) {
            for (int i = 0; i < sizeY; i++) {
                for (int j = 0; j < sizeX; j++) {
                    get(i, j).sub(grid.get(i, j));
                }
            }
            return this;
        }

        public GEMgrid sub(Vec3 v){
            for (int i = 0; i < sizeY; i++) {
                for (int j = 0; j < sizeX; j++) {
                    get(i, j).sub(v);
                }
            }
            return this;
        }

        public GEMgrid normalize(){
            buffer.set(Vec3.Zero);
            final float scale = 1f/grid.length;
            eachTile(v -> buffer.add(v, scale));
            return this;
        }

        public GEMgrid scl(float s) {
            eachTile(v -> v.scl(s));
            return this;
        }

        public GEMgrid addmul(GEMgrid grid, float s){
            for (int i = 0; i < sizeY; i++) {
                for (int j = 0; j < sizeX; j++) {
                    get(i, j).add(grid.get(i, j), s);
                }
            }
            return this;
        }

        public GEMgrid eachTile(Cons<Vec3> cons) {
            for (Vec3 v : grid) {
                cons.get(v);
            }
            return this;
        }
    }
}
