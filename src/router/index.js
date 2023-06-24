import Vue from "vue";
import VueRouter from "vue-router";
import Home from "../views/Home.vue";
import Step from "../views/Step.vue";
import Pipeline from "../views/Pipeline.vue";
import fastqcANDmultiqc from "../views/fastqcANDmultiqc.vue";
import ExpertMode from "../views/ExpertMode.vue";

Vue.use(VueRouter);

const routes = [
  { path: "/step/:stepName/:order", component: Step },
  {
    path: "/home",
    component: Home,
  },
  {
    path: "/premade/:workflowName",
    component: Pipeline,
  },
  {
    path: "/fastqcANDmultiqc",
    component: fastqcANDmultiqc,
  },
  {
    path: "/ExpertMode",
    component: ExpertMode,
  },
];

const router = new VueRouter({
  routes,
});

export default router;
